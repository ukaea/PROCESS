import logging
import math

import numpy as np

from process import process_output as po
from process.exceptions import ProcessValueError
from process.fortran import (
    build_variables,
    buildings_variables,
    constants,
    constraint_variables,
    cost_variables,
    current_drive_variables,
    error_handling,
    fwbs_variables,
    heat_transport_variables,
    numerics,
    pf_power_variables,
    pfcoil_variables,
    physics_variables,
    primary_pumping_variables,
    structure_variables,
    tfcoil_variables,
    times_variables,
)
from process.variables import AnnotatedVariable

logger = logging.getLogger(__name__)


class Power:
    def __init__(self):
        self.outfile = constants.nout

        # Local variables
        self.qmisc = AnnotatedVariable(float, 0.0, docstring="", units="")
        self.qac = AnnotatedVariable(float, 0.0, docstring="", units="")
        self.qcl = AnnotatedVariable(float, 0.0, docstring="", units="")
        self.qss = AnnotatedVariable(float, 0.0, docstring="", units="")
        self.p_shld_coolant_pump_elec_mw = AnnotatedVariable(
            float, 0.0, docstring="", units=""
        )
        self.p_div_coolant_pump_elec_mw = AnnotatedVariable(
            float, 0.0, docstring="", units=""
        )
        self.p_coolant_pump_total_mw = AnnotatedVariable(
            float, 0.0, docstring="", units=""
        )
        self.p_fw_blkt_heat_deposited_mw = AnnotatedVariable(
            float, 0.0, docstring="", units=""
        )
        self.p_fw_blkt_coolant_pump_elec_mw = AnnotatedVariable(
            float, 0.0, docstring="", units=""
        )
        self.p_blkt_breeder_pump_elec_mw = AnnotatedVariable(
            float, 0.0, docstring="", units=""
        )
        self.p_div_heat_deposited_mw = AnnotatedVariable(
            float, 0.0, docstring="", units=""
        )
        self.p_fw_heat_deposited_mw = AnnotatedVariable(
            float, 0.0, docstring="", units=""
        )
        self.p_blkt_heat_deposited_mw = AnnotatedVariable(
            float, 0.0, docstring="", units=""
        )
        self.pthermblkt_liq = AnnotatedVariable(float, 0.0, docstring="", units="")
        self.p_shld_heat_deposited_mw = AnnotatedVariable(
            float, 0.0, docstring="", units=""
        )
        self.p_cp_coolant_pump_elec_mw = AnnotatedVariable(
            float, 0.0, docstring="", units=""
        )
        self.p_plant_core_systems_elec_mw = AnnotatedVariable(
            float, 0.0, docstring="", units=""
        )
        self.pdivfraction = AnnotatedVariable(float, 0.0, docstring="", units="")
        self.delta_eta = AnnotatedVariable(float, 0.0, docstring="", units="")
        self.iprimdiv = AnnotatedVariable(float, 0.0, docstring="", units="")
        self.p_turbine_loss_mw = AnnotatedVariable(float, 0.0, docstring="", units="")

    def pfpwr(self, output: bool):
        """
        PF coil power supply requirements
        author: P J Knight, CCFE, Culham Science Centre
        outfile : input integer : output file unit
        iprint : input integer : switch for writing to output (1=yes)
        This routine calculates the MVA, power and energy requirements
        for the PF coil systems.  Units are MW and MVA for power terms.
        The routine checks at the beginning of the flattop for the
        peak MVA, and at the end of flattop for the peak stored energy.
        The reactive (inductive) components use waves to calculate the
        <I>dI/dt</I> at the time periods.
        None
        """
        powpfii = np.zeros((pfcoil_variables.ngc2,))
        cktr = np.zeros((pfcoil_variables.ngc2,))
        pfcr = np.zeros((pfcoil_variables.ngc2,))
        albusa = np.zeros((pfcoil_variables.ngc2,))
        pfbusr = np.zeros((pfcoil_variables.ngc2,))
        pfcr = np.zeros((pfcoil_variables.ngc2,))
        cktr = np.zeros((pfcoil_variables.ngc2,))
        rcktvm = np.zeros((pfcoil_variables.ngc2,))
        rcktpm = np.zeros((pfcoil_variables.ngc2,))
        vpfi = np.zeros((pfcoil_variables.ngc2,))
        psmva = np.zeros((pfcoil_variables.ngc2,))
        poloidalenergy = np.zeros((6,))
        inductxcurrent = np.zeros((6,))
        pfdissipation = np.zeros((5,))

        #  Bus length
        pfbusl = 8.0e0 * physics_variables.rmajor + 140.0e0

        #  Find power requirements for PF coils at times_variables.tim(ktim)

        #  PF coil resistive power requirements
        #  Bussing losses assume aluminium bussing with 100 A/cm**2
        ic = -1
        ngrpt = pfcoil_variables.n_pf_coil_groups
        if build_variables.iohcl != 0:
            ngrpt = ngrpt + 1

        pf_power_variables.srcktpm = 0.0e0
        pfbuspwr = 0.0e0

        for ig in range(ngrpt):
            ic = ic + pfcoil_variables.n_pf_coils_in_group[ig]

            #  Section area of aluminium bussing for circuit (cm**2)
            #  pfcoil_variables.c_pf_coil_turn_peak_input : max current per turn of coil (A)
            albusa[ig] = abs(pfcoil_variables.c_pf_coil_turn_peak_input[ic]) / 100.0e0

            #  Resistance of bussing for circuit (ohm)
            #  pfbusl : bus length for each PF circuit (m)
            #  pfbusr[ig] = 1.5e0 * 2.62e-4 * pfbusl / albusa[ig]
            #  I have removed the fudge factor of 1.5 but included it in the value of rhopfbus
            pfbusr[ig] = pfcoil_variables.rhopfbus * pfbusl / (albusa[ig] / 10000)

            #  Total PF coil resistance (during burn)
            #  pfcoil_variables.c_pf_cs_coils_peak_ma : maximum current in coil (A)
            pfcr[ig] = (
                pfcoil_variables.rho_pf_coil
                * 2.0e0
                * np.pi
                * pfcoil_variables.r_pf_coil_middle[ic]
                * abs(
                    pfcoil_variables.j_pf_coil_wp_peak[ic]
                    / (
                        (1.0e0 - pfcoil_variables.f_a_pf_coil_void[ic])
                        * 1.0e6
                        * pfcoil_variables.c_pf_cs_coils_peak_ma[ic]
                    )
                )
                * pfcoil_variables.n_pf_coil_turns[ic] ** 2
                * pfcoil_variables.n_pf_coils_in_group[ig]
            )

            cktr[ig] = pfcr[ig] + pfbusr[ig]  # total resistance of circuit (ohms)
            cptburn = (
                pfcoil_variables.c_pf_coil_turn_peak_input[ic]
                * pfcoil_variables.c_pf_cs_coil_pulse_end_ma[ic]
                / pfcoil_variables.c_pf_cs_coils_peak_ma[ic]
            )
            rcktvm[ig] = abs(cptburn) * cktr[ig]  # peak resistive voltage (V)
            rcktpm[ig] = 1.0e-6 * rcktvm[ig] * abs(cptburn)  # peak resistive power (MW)

            #  Compute the sum of resistive power in the PF circuits, kW
            pfbuspwr = pfbuspwr + 1.0e-3 * pfbusr[ig] * cptburn**2
            pf_power_variables.srcktpm = pf_power_variables.srcktpm + 1.0e3 * rcktpm[ig]

        #  Inductive MVA requirements, and stored energy
        delktim = times_variables.t_current_ramp_up

        #  PF system (including Central Solenoid solenoid) inductive MVA requirements
        #  pfcoil_variables.c_pf_coil_turn(i,j) : current per turn of coil i at (end) time period j (A)
        powpfi = 0.0e0
        powpfr = 0.0e0
        powpfr2 = 0.0e0

        #  pfcoil_variables.n_pf_cs_plasma_circuits : total number of PF coils (including Central Solenoid and plasma)
        #          plasma is #n_pf_cs_plasma_circuits, and Central Solenoid is #(pfcoil_variables.n_pf_cs_plasma_circuits-1)
        #  pfcoil_variables.ind_pf_cs_plasma_mutual(i,j) : mutual inductance between coil i and j
        for i in range(pfcoil_variables.n_pf_cs_plasma_circuits):
            powpfii[i] = 0.0e0
            vpfi[i] = 0.0e0

        jpf = -1
        poloidalenergy[:] = 0.0e0
        for jjpf in range(ngrpt):  # Loop over all groups of PF coils.
            for _jjpf2 in range(
                pfcoil_variables.n_pf_coils_in_group[jjpf]
            ):  # Loop over all coils in each group
                jpf = jpf + 1
                inductxcurrent[:] = 0.0e0
                for ipf in range(pfcoil_variables.n_pf_cs_plasma_circuits):
                    #  Voltage in circuit jpf due to change in current from circuit ipf
                    vpfij = (
                        pfcoil_variables.ind_pf_cs_plasma_mutual[jpf, ipf]
                        * (
                            pfcoil_variables.c_pf_coil_turn[ipf, 2]
                            - pfcoil_variables.c_pf_coil_turn[ipf, 1]
                        )
                        / delktim
                    )

                    #  Voltage in circuit jpf at time, times_variables.tim(3), due to changes in coil currents
                    vpfi[jpf] = vpfi[jpf] + vpfij

                    #  MVA in circuit jpf at time, times_variables.tim(3) due to changes in current
                    powpfii[jpf] = (
                        powpfii[jpf]
                        + vpfij * pfcoil_variables.c_pf_coil_turn[jpf, 2] / 1.0e6
                    )

                    # Term used for calculating stored energy at each time
                    for time in range(6):
                        inductxcurrent[time] = (
                            inductxcurrent[time]
                            + pfcoil_variables.ind_pf_cs_plasma_mutual[jpf, ipf]
                            * pfcoil_variables.c_pf_coil_turn[ipf, time]
                        )

                    # engx = engx + pfcoil_variables.ind_pf_cs_plasma_mutual(jpf,ipf)*pfcoil_variables.c_pf_coil_turn(ipf,5)

                #  Stored magnetic energy of the poloidal field at each time
                # 'time' is the time INDEX.  'tim' is the time.
                for time in range(6):
                    poloidalenergy[time] = (
                        poloidalenergy[time]
                        + 0.5e0
                        * inductxcurrent[time]
                        * pfcoil_variables.c_pf_coil_turn[jpf, time]
                    )

                #   do time = 1,5
                #     # Mean rate of change of stored energy between time and time+1
                #     if(abs(times_variables.tim(time+1)-times_variables.tim(time)).gt.1.0e0) :
                #         pf_power_variables.poloidalpower(time) = (poloidalenergy(time+1)-poloidalenergy(time)) / (times_variables.tim(time+1)-times_variables.tim(time))
                #     else:
                #         # Flag when an interval is small or zero MDK 30/11/16
                #         pf_power_variables.poloidalpower(time) = 9.9e9
                #

                #   end do
                #   #engxpc = 0.5e0 * engx * pfcoil_variables.c_pf_coil_turn(jpf,5)
                #   #ensxpf = ensxpf + engxpc

                #  Resistive power in circuits at times times_variables.tim(3) and times_variables.tim(5) respectively (MW)
                powpfr = (
                    powpfr
                    + pfcoil_variables.n_pf_coil_turns[jpf]
                    * pfcoil_variables.c_pf_coil_turn[jpf, 2]
                    * cktr[jjpf]
                    / 1.0e6
                )
                powpfr2 = (
                    powpfr2
                    + pfcoil_variables.n_pf_coil_turns[jpf]
                    * pfcoil_variables.c_pf_coil_turn[jpf, 4]
                    * cktr[jjpf]
                    / 1.0e6
                )
                powpfi = powpfi + powpfii[jpf]

        for time in range(5):
            # Stored magnetic energy of the poloidal field at each time
            # 'time' is the time INDEX.  'tim' is the time.
            # Mean rate of change of stored energy between time and time+1
            if abs(times_variables.tim[time + 1] - times_variables.tim[time]) > 1.0e0:
                pf_power_variables.poloidalpower[time] = (
                    poloidalenergy[time + 1] - poloidalenergy[time]
                ) / (times_variables.tim[time + 1] - times_variables.tim[time])
            else:
                # Flag when an interval is small or zero MDK 30/11/16
                pf_power_variables.poloidalpower[time] = 9.9e9

            # Electrical energy dissipated in PFC power supplies as they increase or decrease the poloidal field energy
            # This assumes that the energy storage in the PFC power supply is lossless and that currents
            # in the coils can be varied without loss when there is no change in the energy in the poloidal field.
            # Energy is dissipated only when energy moves into or out of the store in the power supply.
            # Issue #713
            pfdissipation[time] = abs(
                poloidalenergy[time + 1] - poloidalenergy[time]
            ) * (1.0e0 / pfcoil_variables.etapsu - 1.0e0)

        # Mean power dissipated
        # The flat top duration (time 4 to 5) is the denominator, as this is the time when electricity is generated.
        if times_variables.tim[4] - times_variables.tim[3] > 1.0e0:
            pfpower = sum(pfdissipation[:]) / (
                times_variables.tim[4] - times_variables.tim[3]
            )
        else:
            # Give up when an interval is small or zero.
            pfpower = 0.0e0

        pfpowermw = pfpower / 1.0e6

        #  Compute the maximum stored energy and the maximum dissipative
        #  energy in all the PF circuits over the entire cycle time, MJ
        # ensxpfm = 1.0e-6 * ensxpf
        pf_power_variables.ensxpfm = 1.0e-6 * max(poloidalenergy)
        # Peak absolute rate of change of stored energy in poloidal field (MW)
        pf_power_variables.peakpoloidalpower = (
            max(abs(pf_power_variables.poloidalpower)) / 1.0e6
        )

        #  Maximum total MVA requirements
        heat_transport_variables.peakmva = max((powpfr + powpfi), powpfr2)

        pf_power_variables.vpfskv = 20.0e0
        pf_power_variables.pfckts = (
            pfcoil_variables.n_pf_cs_plasma_circuits - 2
        ) + 6.0e0
        pf_power_variables.spfbusl = pfbusl * pf_power_variables.pfckts
        pf_power_variables.acptmax = 0.0e0
        pf_power_variables.spsmva = 0.0e0

        for jpf in range(pfcoil_variables.n_pf_cs_plasma_circuits - 1):
            #  Power supply MVA for each PF circuit
            psmva[jpf] = 1.0e-6 * abs(
                vpfi[jpf] * pfcoil_variables.c_pf_coil_turn_peak_input[jpf]
            )

            #  Sum of the power supply MVA of the PF circuits
            pf_power_variables.spsmva = pf_power_variables.spsmva + psmva[jpf]

            #  Average of the maximum currents in the PF circuits, kA
            pf_power_variables.acptmax = (
                pf_power_variables.acptmax
                + 1.0e-3
                * abs(pfcoil_variables.c_pf_coil_turn_peak_input[jpf])
                / pf_power_variables.pfckts
            )

        #  PF wall plug power dissipated in power supply for ohmic heating (MW)
        #  This is additional to that required for moving stored energy around
        # pfwpmw = physics_variables.p_plasma_ohmic_mw / pfcoil_variables.etapsu
        wall_plug_ohmicmw = physics_variables.p_plasma_ohmic_mw * (
            1.0e0 / pfcoil_variables.etapsu - 1.0e0
        )
        # Total mean wall plug power dissipated in PFC and CS power supplies.  Issue #713
        pfcoil_variables.pfwpmw = wall_plug_ohmicmw + pfpowermw

        #  Output Section
        if output == 0:
            return

        po.oheadr(self.outfile, "PF Coils and Central Solenoid: Power and Energy")
        po.ovarre(
            self.outfile,
            "Number of PF coil circuits",
            "(pfckts)",
            pf_power_variables.pfckts,
        )
        po.ovarre(
            self.outfile,
            "Sum of PF power supply ratings (MVA)",
            "(spsmva)",
            pf_power_variables.spsmva,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Total PF coil circuit bus length (m)",
            "(spfbusl)",
            pf_power_variables.spfbusl,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Total PF coil bus resistive power (kW)",
            "(pfbuspwr)",
            pfbuspwr,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Total PF coil resistive power (kW)",
            "(srcktpm)",
            pf_power_variables.srcktpm,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Maximum PF coil voltage (kV)",
            "(vpfskv)",
            pf_power_variables.vpfskv,
        )
        po.ovarre(
            self.outfile,
            "Efficiency of transfer of PF stored energy into or out of storage",
            "(etapsu)",
            pfcoil_variables.etapsu,
        )
        po.ocmmnt(
            self.outfile,
            "(Energy is dissipated in PFC power supplies only when total PF energy increases or decreases.)",
        )

        po.ovarre(
            self.outfile,
            "Maximum stored energy in poloidal field (MJ)",
            "(ensxpfm)",
            pf_power_variables.ensxpfm,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Peak absolute rate of change of stored energy in poloidal field (MW)",
            "(peakpoloidalpower)",
            pf_power_variables.peakpoloidalpower,
            "OP ",
        )

        if (numerics.ioptimz > 0) and (numerics.active_constraints[65]):
            po.ovarre(
                self.outfile,
                "Max permitted abs rate of change of stored energy in poloidal field (MW)",
                "maxpoloidalpower",
                pf_power_variables.maxpoloidalpower,
            )

        if any(poloidalenergy < 0.0e0):
            po.oheadr(self.outfile, "ERROR Negative stored energy in poloidal field")
            logger.error(f"{'ERROR Negative stored energy in poloidal field'}")

        po.ocmmnt(self.outfile, "Energy stored in poloidal magnetic field :")
        po.oblnkl(self.outfile)

        # write(self.outfile,50)(times_variables.tim(time),time=1,6)

    def acpow(self, output: bool):
        """
        AC power requirements
        author: P J Knight, CCFE, Culham Science Centre
        author: P C Shipe, ORNL
        outfile : input integer : output file unit
        iprint : input integer : switch for writing to output (1=yes)
        The routine was drastically shortened on 23/01/90 (ORNL) from the
        original TETRA routine to provide only the total power needs for
        the plant. Included in STORAC in January 1992 by P.C. Shipe.
        None
        """
        ptfmw = heat_transport_variables.tfacpd

        # Power to PF coil power supplies, MW
        ppfmw = 1.0e-3 * pf_power_variables.srcktpm

        if pf_power_variables.iscenr == 2:
            ppfmw = ppfmw + heat_transport_variables.peakmva

        #  Power to plasma heating supplies, MW
        pheatingmw = (
            heat_transport_variables.p_hcd_electric_total_mw
        )  # Should be zero if i_plasma_ignited==1

        #  Power to cryogenic comp. motors, MW
        crymw = heat_transport_variables.p_cryo_plant_electric_mw

        #  Facility base load, MW (loads not dependent on floor area)
        basemw = heat_transport_variables.p_plant_electric_base * 1.0e-6

        #  Power needed per unit floor area, kW/m2
        pkwpm2 = heat_transport_variables.pwpm2 * 1.0e-3

        #  Power to divertor coil supplies, MW
        bdvmw = 0.0e0

        #  Total pulsed power system load, MW
        heat_transport_variables.pacpmw = (
            ppfmw
            + bdvmw
            + ptfmw
            + crymw
            + heat_transport_variables.vachtmw
            + heat_transport_variables.p_coolant_pump_elec_total_mw
            + heat_transport_variables.p_tritium_plant_electric_mw
            + pheatingmw
        )

        #  Add contribution from motor-generator flywheels if these are part of
        #  the PF coil energy storage system
        if pf_power_variables.iscenr != 2:
            heat_transport_variables.pacpmw = (
                heat_transport_variables.pacpmw + heat_transport_variables.fmgdmw
            )

        #  Total baseline power to facility loads, MW
        heat_transport_variables.fcsht = (
            basemw + buildings_variables.efloor * pkwpm2 / 1000.0e0
        )

        # Estimate of the total low voltage power, MW
        # MDK No idea what this is - especially the last term
        # It is used in the old cost routine, so I will leave it in place.
        heat_transport_variables.tlvpmw = (
            heat_transport_variables.fcsht
            + heat_transport_variables.p_tritium_plant_electric_mw
            + heat_transport_variables.p_coolant_pump_elec_total_mw
            + heat_transport_variables.vachtmw
            + 0.5e0 * (crymw + ppfmw)
        )

        if output == 0:
            return

        #  Output section
        # po.oheadr(self.outfile,'AC Power')
        po.oheadr(self.outfile, "Electric Power Requirements")
        po.ovarre(self.outfile, "Facility base load (MW)", "(basemw)", basemw)
        po.ovarre(self.outfile, "Divertor coil power supplies (MW)", "(bdvmw)", bdvmw)
        po.ovarre(
            self.outfile, "Cryoplant electric power (MW)", "(crymw)", crymw, "OP "
        )
        # po.ovarre(self.outfile,'Heat removed from cryogenic coils (MWth)','(helpow/1.0e6)',helpow/1.0e6)
        # po.ovarre(self.outfile,'MGF (motor-generator flywheel) units (MW)', '(fmgdmw)',fmgdmw)
        # po.ovarin(self.outfile,'Primary coolant pumps (MW)', '(i_blkt_coolant_type)',i_blkt_coolant_type)
        po.ovarre(
            self.outfile,
            "Primary coolant pumps (MW)",
            "(p_coolant_pump_elec_total_mw..)",
            heat_transport_variables.p_coolant_pump_elec_total_mw,
            "OP ",
        )

        po.ovarre(self.outfile, "PF coil power supplies (MW)", "(ppfmw)", ppfmw, "OP ")
        # po.ovarre(self.outfile,'Power/floor area (kW/m2)','(pkwpm2)',pkwpm2)
        po.ovarre(self.outfile, "TF coil power supplies (MW)", "(ptfmw)", ptfmw, "OP ")
        po.ovarre(
            self.outfile,
            "Plasma heating supplies (MW)",
            "(pheatingmw)",
            pheatingmw,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Tritium processing (MW)",
            "(p_tritium_plant_electric_mw..)",
            heat_transport_variables.p_tritium_plant_electric_mw,
        )
        po.ovarre(
            self.outfile,
            "Vacuum pumps  (MW)",
            "(vachtmw..)",
            heat_transport_variables.vachtmw,
        )

        po.oblnkl(self.outfile)

        po.ovarre(
            self.outfile,
            "Total pulsed power (MW)",
            "(pacpmw)",
            heat_transport_variables.pacpmw,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Total base power required at all times (MW)",
            "(fcsht)",
            heat_transport_variables.fcsht,
            "OP ",
        )
        # MDK Remove this output: no idea what this is
        # po.ovarre(self.outfile,'Total low voltage power (MW)','(tlvpmw)',tlvpmw)

    def power1(self):
        """
        Calculates the first part of the heat transport
        and plant power balance constituents
        author: P J Knight, CCFE, Culham Science Centre
        None
        This routine calculates the first part of the heat transport
        and plant power balance constituents.
        None
        """
        if (
            fwbs_variables.i_coolant_pumping != 2
            and fwbs_variables.i_coolant_pumping != 3
        ):
            primary_pumping_variables.p_fw_blkt_coolant_pump_mw = (
                heat_transport_variables.p_fw_coolant_pump_mw
                + heat_transport_variables.p_blkt_coolant_pump_mw
            )

        #  Account for pump electrical inefficiencies. The coolant pumps are not assumed to be
        #  100% efficient so the electric power to run them is greater than the power deposited
        #  in the coolant.  The difference should be lost as secondary heat.
        self.p_fw_blkt_coolant_pump_elec_mw = (
            primary_pumping_variables.p_fw_blkt_coolant_pump_mw / fwbs_variables.etahtp
        )
        self.p_shld_coolant_pump_elec_mw = (
            heat_transport_variables.p_shld_coolant_pump_mw / fwbs_variables.etahtp
        )
        self.p_div_coolant_pump_elec_mw = (
            heat_transport_variables.p_div_coolant_pump_mw / fwbs_variables.etahtp
        )
        if (
            fwbs_variables.i_blkt_dual_coolant > 0
            and fwbs_variables.i_coolant_pumping == 2
        ):
            self.p_blkt_breeder_pump_elec_mw = (
                heat_transport_variables.p_blkt_breeder_pump_mw / fwbs_variables.etahtp
            )

        if (
            fwbs_variables.i_blkt_dual_coolant > 0
            and fwbs_variables.i_coolant_pumping == 2
        ):
            # Total mechanical pump power (deposited in coolant)
            self.p_coolant_pump_total_mw = (
                primary_pumping_variables.p_fw_blkt_coolant_pump_mw
                + heat_transport_variables.p_blkt_breeder_pump_mw
                + heat_transport_variables.p_shld_coolant_pump_mw
                + heat_transport_variables.p_div_coolant_pump_mw
            )
            # Minimum total electrical power for primary coolant pumps  (MW) Issue #303
            # Recommended to leave the minimum value at zero.
            # Note that p_coolant_pump_elec_total_mw is an ELECTRICAL power
            heat_transport_variables.p_coolant_pump_elec_total_mw = (
                self.p_fw_blkt_coolant_pump_elec_mw
                + self.p_blkt_breeder_pump_elec_mw
                + self.p_shld_coolant_pump_elec_mw
                + self.p_div_coolant_pump_elec_mw
            )

        else:
            # Total mechanical pump power (deposited in coolant)
            self.p_coolant_pump_total_mw = (
                primary_pumping_variables.p_fw_blkt_coolant_pump_mw
                + heat_transport_variables.p_shld_coolant_pump_mw
                + heat_transport_variables.p_div_coolant_pump_mw
            )

            # Minimum total electrical power for primary coolant pumps  (MW) Issue #303
            # Recommended to leave the minimum value at zero.
            # Note that p_coolant_pump_elec_total_mw is an ELECTRICAL power
            heat_transport_variables.p_coolant_pump_elec_total_mw = (
                self.p_fw_blkt_coolant_pump_elec_mw
                + self.p_shld_coolant_pump_elec_mw
                + self.p_div_coolant_pump_elec_mw,
            )

        #  Heat lost through pump power inefficiencies (MW)
        heat_transport_variables.p_coolant_pump_loss_total_mw = (
            heat_transport_variables.p_coolant_pump_elec_total_mw
            - self.p_coolant_pump_total_mw
        )

        # Calculate total deposited power (MW), n.b. energy multiplication in p_blkt_nuclear_heat_total_mw already

        if fwbs_variables.i_coolant_pumping == 2:
            # Liquid metal breeder/coolant
            if fwbs_variables.i_blkt_dual_coolant == 2:
                self.pthermblkt_liq = (
                    fwbs_variables.p_blkt_nuclear_heat_total_mw
                    * fwbs_variables.f_nuc_pow_bz_liq
                ) + heat_transport_variables.p_blkt_breeder_pump_mw
            elif fwbs_variables.i_blkt_dual_coolant == 1:
                self.pthermblkt_liq = heat_transport_variables.p_blkt_breeder_pump_mw

            # First wall and blanket coolant combined
            if fwbs_variables.i_blkt_dual_coolant == 2:
                self.p_fw_blkt_heat_deposited_mw = (
                    self.pthermblkt_liq
                    + fwbs_variables.p_fw_nuclear_heat_total_mw
                    + fwbs_variables.p_fw_rad_total_mw
                    + (
                        fwbs_variables.p_blkt_nuclear_heat_total_mw
                        * (1 - fwbs_variables.f_nuc_pow_bz_liq)
                    )
                    + primary_pumping_variables.p_fw_blkt_coolant_pump_mw
                    + current_drive_variables.p_beam_orbit_loss_mw
                    + physics_variables.p_fw_alpha_mw
                    + current_drive_variables.p_beam_shine_through_mw
                )
            elif fwbs_variables.i_blkt_dual_coolant == 1:
                self.p_fw_blkt_heat_deposited_mw = (
                    self.pthermblkt_liq
                    + fwbs_variables.p_fw_nuclear_heat_total_mw
                    + fwbs_variables.p_fw_rad_total_mw
                    + fwbs_variables.p_blkt_nuclear_heat_total_mw
                    + primary_pumping_variables.p_fw_blkt_coolant_pump_mw
                    + current_drive_variables.p_beam_orbit_loss_mw
                    + physics_variables.p_fw_alpha_mw
                    + current_drive_variables.p_beam_shine_through_mw
                )
            else:
                self.p_fw_blkt_heat_deposited_mw = (
                    fwbs_variables.p_fw_nuclear_heat_total_mw
                    + fwbs_variables.p_fw_rad_total_mw
                    + fwbs_variables.p_blkt_nuclear_heat_total_mw
                    + primary_pumping_variables.p_fw_blkt_coolant_pump_mw
                    + current_drive_variables.p_beam_orbit_loss_mw
                    + physics_variables.p_fw_alpha_mw
                    + current_drive_variables.p_beam_shine_through_mw
                )

        elif fwbs_variables.i_coolant_pumping == 3:
            # First wall and blanket coolant combined
            self.p_fw_blkt_heat_deposited_mw = (
                fwbs_variables.p_fw_nuclear_heat_total_mw
                + fwbs_variables.p_fw_rad_total_mw
                + fwbs_variables.p_blkt_nuclear_heat_total_mw
                + primary_pumping_variables.p_fw_blkt_coolant_pump_mw
                + current_drive_variables.p_beam_orbit_loss_mw
                + physics_variables.p_fw_alpha_mw
                + current_drive_variables.p_beam_shine_through_mw
            )

        else:
            #  Total power deposited in first wall coolant (MW)
            self.p_fw_heat_deposited_mw = (
                fwbs_variables.p_fw_nuclear_heat_total_mw
                + fwbs_variables.p_fw_rad_total_mw
                + heat_transport_variables.p_fw_coolant_pump_mw
                + current_drive_variables.p_beam_orbit_loss_mw
                + physics_variables.p_fw_alpha_mw
                + current_drive_variables.p_beam_shine_through_mw
            )
            #  Total power deposited in blanket coolant (MW) (energy multiplication in fwbs_variables.p_blkt_nuclear_heat_total_mw already)
            self.p_blkt_heat_deposited_mw = (
                fwbs_variables.p_blkt_nuclear_heat_total_mw
                + heat_transport_variables.p_blkt_coolant_pump_mw
            )
            self.p_fw_blkt_heat_deposited_mw = (
                self.p_fw_heat_deposited_mw + self.p_blkt_heat_deposited_mw
            )

        #  Total power deposited in shield coolant (MW)
        self.p_shld_heat_deposited_mw = (
            fwbs_variables.p_cp_shield_nuclear_heat_mw
            + fwbs_variables.p_shld_nuclear_heat_mw
            + heat_transport_variables.p_shld_coolant_pump_mw
        )

        #  Total thermal power deposited in divertor coolant (MW)
        #  = (conduction to divertor, less radiation) + (neutron and radiation power)
        #  using physics_variables.p_plasma_separatrix_mw as calculated in physics.f90
        self.p_div_heat_deposited_mw = (
            physics_variables.p_plasma_separatrix_mw
            + (
                fwbs_variables.p_div_nuclear_heat_total_mw
                + fwbs_variables.p_div_rad_total_mw
            )
            + heat_transport_variables.p_div_coolant_pump_mw
        )

        #  Heat removal from first wall and divertor (MW) (only used in costs.f90)
        if fwbs_variables.i_coolant_pumping != 3:
            heat_transport_variables.p_fw_div_heat_deposited_mw = (
                self.p_fw_heat_deposited_mw + self.p_div_heat_deposited_mw
            )

        #  Thermal to electric efficiency
        heat_transport_variables.eta_turbine = self.plant_thermal_efficiency(
            heat_transport_variables.eta_turbine
        )
        heat_transport_variables.etath_liq = self.plant_thermal_efficiency_2(
            heat_transport_variables.etath_liq
        )

        #  Primary (high-grade) thermal power, available for electricity generation.  Switch heat_transport_variables.iprimshld
        #  is 1 or 0, is user choice on whether the shield thermal power goes to primary or secondary heat
        if fwbs_variables.i_thermal_electric_conversion == 0:
            #  Primary thermal power (MW)
            heat_transport_variables.pthermmw = (
                self.p_fw_blkt_heat_deposited_mw
                + heat_transport_variables.iprimshld * self.p_shld_heat_deposited_mw
            )
            #  Secondary thermal power deposited in divertor (MW)
            heat_transport_variables.psecdiv = self.p_div_heat_deposited_mw
            # Divertor primary/secondary power switch: does NOT contribute to energy generation cycle
            self.iprimdiv = 0
        else:
            #  Primary thermal power (MW)
            heat_transport_variables.pthermmw = (
                self.p_fw_blkt_heat_deposited_mw
                + heat_transport_variables.iprimshld * self.p_shld_heat_deposited_mw
                + self.p_div_heat_deposited_mw
            )
            #  Secondary thermal power deposited in divertor (MW)
            heat_transport_variables.psecdiv = 0.0e0
            # Divertor primary/secondary power switch: contributes to energy generation cycle
            self.iprimdiv = 1

        if abs(heat_transport_variables.pthermmw) < 1.0e-4:
            logger.error(f"{'ERROR Primary thermal power is zero or negative'}")

        # #284 Fraction of total high-grade thermal power to divertor
        self.pdivfraction = (
            self.p_div_heat_deposited_mw / heat_transport_variables.pthermmw
        )
        # Loss in efficiency as this primary power is collecetd at very low temperature
        self.delta_eta = 0.339 * self.pdivfraction

        #  Secondary thermal power deposited in shield
        heat_transport_variables.psecshld = self.p_shld_heat_deposited_mw * (
            1 - heat_transport_variables.iprimshld
        )

        #  Secondary thermal power lost to HCD apparatus and diagnostics
        heat_transport_variables.psechcd = (
            fwbs_variables.p_fw_hcd_nuclear_heat_mw
            + fwbs_variables.p_fw_hcd_rad_total_mw
        )

        #  Number of primary heat exchangers
        heat_transport_variables.nphx = math.ceil(
            heat_transport_variables.pthermmw / 1000.0e0
        )

        #  Secondary heat (some of it... rest calculated in POWER2)
        #  Wall plug injection power
        # MDK
        # heat_transport_variables.p_hcd_electric_total_mw = (current_drive_variables.p_hcd_injected_total_mw + current_drive_variables.p_beam_orbit_loss_mw + physics_variables.p_fw_alpha_mw)/eta_hcd_primary_injector_wall_plug
        # heat_transport_variables.p_hcd_electric_total_mw calculated in current_drive.f90

        #  Waste injection power
        if physics_variables.i_plasma_ignited == 0:
            # MDK
            # pinjht = heat_transport_variables.p_hcd_electric_total_mw - current_drive_variables.p_hcd_injected_total_mw - current_drive_variables.p_beam_orbit_loss_mw - physics_variables.p_fw_alpha_mw
            heat_transport_variables.pinjht = (
                heat_transport_variables.p_hcd_electric_total_mw
                - current_drive_variables.p_hcd_injected_total_mw
            )
        else:
            heat_transport_variables.pinjht = 0.0e0

        #  Cryogenic power
        # ---
        # Initialisation (unchanged if all coil resisitive)
        heat_transport_variables.helpow = 0.0e0
        heat_transport_variables.p_cryo_plant_electric_mw = 0.0e0
        p_tf_cryoal_cryo = 0.0e0
        tfcoil_variables.cryo_cool_req = 0.0e0

        # Superconductors TF/PF cryogenic cooling
        if tfcoil_variables.i_tf_sup == 1 or pfcoil_variables.i_pf_conductor == 0:
            # heat_transport_variables.helpow calculation
            heat_transport_variables.helpow = self.cryo(
                tfcoil_variables.i_tf_sup,
                tfcoil_variables.tfcryoarea,
                structure_variables.coldmass,
                fwbs_variables.ptfnuc,
                pf_power_variables.ensxpfm,
                times_variables.t_pulse_repetition,
                tfcoil_variables.cpttf,
                tfcoil_variables.n_tf_coils,
            )

            # Use 13% of ideal Carnot efficiency to fit J. Miller estimate
            # Rem SK : This ITER efficiency is very low compare to the Strowbridge curve
            #          any reasons why?
            # Calculate electric power requirement for cryogenic plant at tfcoil_variables.temp_tf_cryo (MW)
            heat_transport_variables.p_cryo_plant_electric_mw = (
                1.0e-6
                * (constants.temp_room - tfcoil_variables.temp_tf_cryo)
                / (tfcoil_variables.eff_tf_cryo * tfcoil_variables.temp_tf_cryo)
                * heat_transport_variables.helpow
            )

        # Cryogenic alumimium
        # Rem : The carnot efficiency is assumed at 40% as this is a conservative assumption since a 50%
        #       has been deduced from detailed studies
        # Rem : Nuclear heating on the outer legs assumed to be negligible
        # Rem : To be updated with 2 cooling loops for TART designs
        if tfcoil_variables.i_tf_sup == 2:
            # Heat removal power at cryogenic temperature tfcoil_variables.tcoolin (W)
            heat_transport_variables.helpow_cryal = (
                tfcoil_variables.p_cp_resistive
                + tfcoil_variables.p_tf_leg_resistive
                + tfcoil_variables.pres_joints
                + fwbs_variables.pnuc_cp_tf * 1.0e6
            )

            # Calculate electric power requirement for cryogenic plant at tfcoil_variables.tcoolin (MW)
            p_tf_cryoal_cryo = (
                1.0e-6
                * (constants.temp_room - tfcoil_variables.tcoolin)
                / (tfcoil_variables.eff_tf_cryo * tfcoil_variables.tcoolin)
                * heat_transport_variables.helpow_cryal
            )

            # Add to electric power requirement for cryogenic plant (MW)
            heat_transport_variables.p_cryo_plant_electric_mw = (
                heat_transport_variables.p_cryo_plant_electric_mw + p_tf_cryoal_cryo
            )

        # Calculate cryo cooling requirement at 4.5K (kW)
        tfcoil_variables.cryo_cool_req = (
            heat_transport_variables.helpow
            * ((293 / tfcoil_variables.temp_tf_cryo) - 1)
            / ((293 / 4.5) - 1)
            + heat_transport_variables.helpow_cryal
            * ((293 / tfcoil_variables.tcoolin) - 1)
            / ((293 / 4.5) - 1)
        ) / 1.0e3

    def power2(self, output: bool):
        """
        Calculates the remainder of the heat transport
        and plant power balance constituents
        author: P J Knight, CCFE, Culham Science Centre
        outfile : input integer : output file unit
        iprint : input integer : switch for writing to output (1=yes)
        This routine calculates the rest of the heat transport
        and plant power balance constituents, not already calculated in
        <A HREF="acpow.html">ACPOW</A> or <A HREF="power1.html">POWER1</A>.
        None
        """
        if physics_variables.itart == 1 and tfcoil_variables.i_tf_sup == 0:
            self.p_cp_coolant_pump_elec_mw = (
                1.0e-6 * tfcoil_variables.p_cp_coolant_pump_elec
            )
        else:
            self.p_cp_coolant_pump_elec_mw = 0.0e0

        #  Facility heat removal (heat_transport_variables.fcsht calculated in ACPOW)
        heat_transport_variables.fachtmw = heat_transport_variables.fcsht

        #  Electrical power consumed by fusion power core systems
        #  (excluding heat transport pumps and auxiliary injection power system)
        #  pfcoil_variables.pfwpmw = Mean electrical energy dissipated in PFC power supplies as they
        #  increase or decrease the poloidal field energy AND extra due to ohmic heating
        #  of the plasma.  Issue #713
        self.p_plant_core_systems_elec_mw = (
            heat_transport_variables.p_cryo_plant_electric_mw
            + heat_transport_variables.fachtmw
            + self.p_cp_coolant_pump_elec_mw
            + heat_transport_variables.tfacpd
            + heat_transport_variables.p_tritium_plant_electric_mw
            + heat_transport_variables.vachtmw
            + pfcoil_variables.pfwpmw
        )

        #  Total secondary heat
        #  (total low-grade heat rejected - does not contribute to power conversion cycle)
        #  Included fwbs_variables.ptfnuc
        # p_plant_secondary_heat_mw = self.p_plant_core_systems_elec_mw + heat_transport_variables.pinjht + heat_transport_variables.p_coolant_pump_loss_total_mw + hthermmw + heat_transport_variables.psecdiv + heat_transport_variables.psecshld + heat_transport_variables.psechcd + fwbs_variables.ptfnuc
        heat_transport_variables.p_plant_secondary_heat_mw = (
            self.p_plant_core_systems_elec_mw
            + heat_transport_variables.pinjht
            + heat_transport_variables.p_coolant_pump_loss_total_mw
            + heat_transport_variables.psecdiv
            + heat_transport_variables.psecshld
            + heat_transport_variables.psechcd
            + fwbs_variables.ptfnuc
        )

        #  Calculate powers relevant to a power-producing plant
        if cost_variables.ireactor == 1:
            #  Gross electric power
            # p_plant_electric_gross_mw = (heat_transport_variables.pthermmw-hthermmw) * heat_transport_variables.eta_turbine
            if (
                fwbs_variables.i_blkt_dual_coolant > 0
                and fwbs_variables.i_coolant_pumping == 2
            ):
                heat_transport_variables.p_plant_electric_gross_mw = (
                    (heat_transport_variables.pthermmw - self.pthermblkt_liq)
                    * heat_transport_variables.eta_turbine
                    + self.pthermblkt_liq * heat_transport_variables.etath_liq
                )
            else:
                heat_transport_variables.p_plant_electric_gross_mw = (
                    heat_transport_variables.pthermmw
                    * heat_transport_variables.eta_turbine
                )

            #  Total recirculating power
            heat_transport_variables.p_plant_electric_recirc_mw = (
                self.p_plant_core_systems_elec_mw
                + heat_transport_variables.p_hcd_electric_total_mw
                + heat_transport_variables.p_coolant_pump_elec_total_mw
            )

            #  Net electric power
            heat_transport_variables.p_plant_electric_net_mw = (
                heat_transport_variables.p_plant_electric_gross_mw
                - heat_transport_variables.p_plant_electric_recirc_mw
            )

            #  Recirculating power fraction
            cirpowfr = (
                heat_transport_variables.p_plant_electric_gross_mw
                - heat_transport_variables.p_plant_electric_net_mw
            ) / heat_transport_variables.p_plant_electric_gross_mw

        if output == 0:
            return

        # TODO: Can output unphysical values if there are no cryogenics - could be omitted from OUT.DAT in this case but leave in for MFILE?
        #  Output section
        po.oheadr(self.outfile, "Cryogenics")
        po.ovarre(
            self.outfile,
            "Conduction and radiation heat loads on cryogenic components (MW)",
            "(qss/1.0d6)",
            self.qss / 1.0e6,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Nuclear heating of cryogenic components (MW)",
            "(qnuc/1.0d6)",
            fwbs_variables.qnuc / 1.0e6,
            "OP ",
        )
        if fwbs_variables.inuclear == 1:
            po.ocmmnt(
                self.outfile, "Nuclear heating of cryogenic components is a user input."
            )
        po.ovarre(
            self.outfile,
            "AC losses in cryogenic components (MW)",
            "(qac/1.0d6)",
            self.qac / 1.0e6,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Resistive losses in current leads (MW)",
            "(qcl/1.0d6)",
            self.qcl / 1.0e6,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "45% allowance for heat loads in transfer lines, storage tanks etc (MW)",
            "(qmisc/1.0d6)",
            self.qmisc / 1.0e6,
            "OP ",
        )

        po.ovarre(
            self.outfile,
            "Sum = Total heat removal at cryogenic temperatures (temp_tf_cryo & tcoolin) (MW)",
            "(helpow + helpow_cryal/1.0d6)",
            (heat_transport_variables.helpow + heat_transport_variables.helpow_cryal)
            * 1.0e-6,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Temperature of cryogenic superconducting components (K)",
            "(temp_tf_cryo)",
            tfcoil_variables.temp_tf_cryo,
        )
        po.ovarre(
            self.outfile,
            "Temperature of cryogenic aluminium components (K)",
            "(tcoolin)",
            tfcoil_variables.tcoolin,
        )
        # TODO: Both of these efficiencies are printed when it should be either 13% (ITER) or 40% (Strawbrige) - subset of TODO on line 1118
        po.ovarre(
            self.outfile,
            "Efficiency (figure of merit) of cryogenic plant is 13% of ideal Carnot value:",
            "",
            (tfcoil_variables.eff_tf_cryo * tfcoil_variables.temp_tf_cryo)
            / (constants.temp_room - tfcoil_variables.temp_tf_cryo),
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Efficiency (figure of merit) of cryogenic aluminium plant is 40% of ideal Carnot value:",
            "",
            (tfcoil_variables.eff_tf_cryo * tfcoil_variables.tcoolin)
            / (constants.temp_room - tfcoil_variables.tcoolin),
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Electric power for cryogenic plant (MW)",
            "(p_cryo_plant_electric_mw)",
            heat_transport_variables.p_cryo_plant_electric_mw,
            "OP ",
        )

        po.oheadr(self.outfile, "Plant Power / Heat Transport Balance")
        if heat_transport_variables.p_plant_electric_net_mw < 0:
            po.ocmmnt(
                self.outfile, "WARNING: Calculated net electric power is negative"
            )
            po.ocmmnt(
                self.outfile, "--------------------------------------------------"
            )

        po.osubhd(self.outfile, "Assumptions :")

        po.ovarre(
            self.outfile,
            "Neutron power multiplication in blanket",
            "(emult)",
            fwbs_variables.emult,
        )

        if physics_variables.n_divertors == 2:
            # Double null configuration
            po.ovarre(
                self.outfile,
                "Double Null Divertor area fraction of whole toroid surface",
                "(2*f_ster_div_single)",
                2.0e0 * fwbs_variables.f_ster_div_single,
            )
        else:
            # Single null configuration
            po.ovarre(
                self.outfile,
                "Divertor area fraction of whole toroid surface",
                "(f_ster_div_single)",
                fwbs_variables.f_ster_div_single,
            )

        po.ovarre(
            self.outfile,
            "H/CD apparatus + diagnostics area fraction",
            "(f_a_fw_hcd)",
            fwbs_variables.f_a_fw_hcd,
        )

        if physics_variables.n_divertors == 2:
            # Double null configuration
            po.ovarre(
                self.outfile,
                "First wall area fraction ",
                "(1-2fdiv-f_a_fw_hcd)",
                1.0e0
                - 2.0e0 * fwbs_variables.f_ster_div_single
                - fwbs_variables.f_a_fw_hcd,
            )
        else:
            # Single null configuration
            po.ovarre(
                self.outfile,
                "First wall area fraction ",
                "(1-f_ster_div_single-f_a_fw_hcd)",
                1.0e0 - fwbs_variables.f_ster_div_single - fwbs_variables.f_a_fw_hcd,
            )

        po.ovarin(
            self.outfile,
            "Switch for pumping of primary coolant",
            "(i_coolant_pumping)",
            fwbs_variables.i_coolant_pumping,
        )
        if fwbs_variables.i_coolant_pumping == 0:
            po.ocmmnt(self.outfile, "User sets mechanical pumping power directly")
        elif fwbs_variables.i_coolant_pumping == 1:
            po.ocmmnt(
                self.outfile,
                "User sets mechanical pumping power as a fraction of thermal power removed by coolant",
            )
        elif fwbs_variables.i_coolant_pumping == 2:
            po.ocmmnt(
                self.outfile,
                "Mechanical pumping power is calculated for FW and blanket",
            )
        elif fwbs_variables.i_coolant_pumping == 3:
            po.ocmmnt(
                self.outfile, "Mechanical pumping power for FW and blanket cooling loop"
            )
            po.ocmmnt(
                self.outfile, "includes heat exchanger, using specified pressure drop"
            )

        po.ovarre(
            self.outfile,
            "Mechanical pumping power for FW cooling loop including heat exchanger (MW)",
            "(p_fw_coolant_pump_mw)",
            heat_transport_variables.p_fw_coolant_pump_mw,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Mechanical pumping power for blanket cooling loop including heat exchanger (MW)",
            "(p_blkt_coolant_pump_mw)",
            heat_transport_variables.p_blkt_coolant_pump_mw,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Mechanical pumping power for FW and blanket cooling loop including heat exchanger (MW)",
            "(p_fw_blkt_coolant_pump_mw)",
            primary_pumping_variables.p_fw_blkt_coolant_pump_mw,
            "OP ",
        )

        if fwbs_variables.i_coolant_pumping != 3:
            po.ovarre(
                self.outfile,
                "Mechanical pumping power for FW (MW)",
                "(p_fw_coolant_pump_mw)",
                heat_transport_variables.p_fw_coolant_pump_mw,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Mechanical pumping power for blanket (MW)",
                "(p_blkt_coolant_pump_mw)",
                heat_transport_variables.p_blkt_coolant_pump_mw,
                "OP ",
            )

        po.ovarre(
            self.outfile,
            "Mechanical pumping power for divertor (MW)",
            "(p_div_coolant_pump_mw)",
            heat_transport_variables.p_div_coolant_pump_mw,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Mechanical pumping power for shield and vacuum vessel (MW)",
            "(p_shld_coolant_pump_mw)",
            heat_transport_variables.p_shld_coolant_pump_mw,
            "OP ",
        )

        po.ovarre(
            self.outfile,
            "Electrical pumping power for FW and blanket (MW)",
            "(p_fw_blkt_coolant_pump_elec_mw)",
            self.p_fw_blkt_coolant_pump_elec_mw,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Electrical pumping power for shield (MW)",
            "(p_shld_coolant_pump_elec_mw)",
            self.p_shld_coolant_pump_elec_mw,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Electrical pumping power for divertor (MW)",
            "(p_div_coolant_pump_elec_mw)",
            self.p_div_coolant_pump_elec_mw,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Total electrical pumping power for primary coolant (MW)",
            "(p_coolant_pump_elec_total_mw)",
            heat_transport_variables.p_coolant_pump_elec_total_mw,
            "OP ",
        )

        if fwbs_variables.i_coolant_pumping == 1:
            po.ovarre(
                self.outfile,
                "Coolant pump power / non-pumping thermal power in first wall",
                "(fpumpfw)",
                heat_transport_variables.fpumpfw,
            )
            po.ovarre(
                self.outfile,
                "Coolant pump power / non-pumping thermal power in blanket",
                "(fpumpblkt)",
                heat_transport_variables.fpumpblkt,
            )

        if fwbs_variables.i_coolant_pumping != 0:
            po.ovarre(
                self.outfile,
                "Coolant pump power / non-pumping thermal power in shield",
                "(fpumpshld)",
                heat_transport_variables.fpumpshld,
            )
            po.ovarre(
                self.outfile,
                "Coolant pump power / non-pumping thermal power in divertor",
                "(fpumpdiv)",
                heat_transport_variables.fpumpdiv,
            )

        po.ovarre(
            self.outfile,
            "Electrical efficiency of heat transport coolant pumps",
            "(etahtp)",
            fwbs_variables.etahtp,
        )
        # #284
        po.osubhd(self.outfile, "Plant thermodynamics: options :")
        if self.iprimdiv == 1:
            po.ocmmnt(
                self.outfile,
                "Divertor thermal power is collected at only 150 C and is used to preheat the coolant in the power cycle",
            )
        elif self.iprimdiv == 0:
            po.ocmmnt(
                self.outfile,
                "Divertor thermal power is not used, but rejected directly to the environment.",
            )

        if heat_transport_variables.iprimshld == 1:
            po.ocmmnt(
                self.outfile,
                "Shield thermal power is collected at only 150 C and is used to preheat the coolant in the power cycle",
            )
        elif heat_transport_variables.iprimshld == 0:
            po.ocmmnt(
                self.outfile,
                "Shield thermal power is not used, but rejected directly to the environment.",
            )

        if cost_variables.ireactor == 1:
            if fwbs_variables.i_thermal_electric_conversion == 0:
                po.ocmmnt(
                    self.outfile,
                    "Power conversion cycle efficiency model: "
                    "efficiency set according to blanket type (div power to secondary)",
                )
            elif fwbs_variables.i_thermal_electric_conversion == 1:
                po.ocmmnt(
                    self.outfile,
                    "Power conversion cycle efficiency model: "
                    "efficiency set according to blanket type (div power to primary)",
                )
                po.ovarrf(
                    self.outfile,
                    "Thermal to electric conversion efficiency of the power conversion cycle",
                    "(eta_turbine)",
                    heat_transport_variables.eta_turbine,
                )
            elif fwbs_variables.i_thermal_electric_conversion == 2:
                po.ocmmnt(
                    self.outfile,
                    "Power conversion cycle efficiency model: user-defined efficiency",
                )
                po.ovarrf(
                    self.outfile,
                    "Thermal to electric conversion efficiency of the power conversion cycle",
                    "(eta_turbine)",
                    heat_transport_variables.eta_turbine,
                )
            elif fwbs_variables.i_thermal_electric_conversion == 3:
                po.ocmmnt(
                    self.outfile,
                    "Power conversion cycle efficiency model: steam Rankine cycle",
                )
            else:
                po.ocmmnt(
                    self.outfile,
                    "Power conversion cycle efficiency model: supercritical CO2 cycle",
                )

            if fwbs_variables.i_thermal_electric_conversion > 2:
                po.ovarrf(
                    self.outfile,
                    "Coolant temperature at turbine inlet (K)",
                    "(temp_turbine_coolant_in)",
                    heat_transport_variables.temp_turbine_coolant_in,
                )

            po.ovarrf(
                self.outfile,
                "Fraction of total high-grade thermal power to divertor",
                "(pdivfraction)",
                self.pdivfraction,
                "OP ",
            )

        po.oblnkl(self.outfile)
        po.ocmmnt(
            self.outfile,
            "Power Balance for Reactor (across vacuum vessel boundary) - Detail",
        )
        po.ocmmnt(
            self.outfile,
            "------------------------------------------------------------------",
        )

        pinj = (
            current_drive_variables.p_hcd_injected_total_mw
            if physics_variables.i_plasma_ignited == 0
            else 0.0
        )

        primsum = 0.0e0
        secsum = 0.0e0

        po.oblnkl(self.outfile)
        po.write(self.outfile, "High-grade             Low-grade              Total")
        po.write(self.outfile, "thermal power (MW)     thermal power (MW)      (MW)")

        po.write(self.outfile, "First wall:")
        po.dblcol(
            self.outfile,
            "p_fw_nuclear_heat_total_mw",
            0.0e0,
            fwbs_variables.p_fw_nuclear_heat_total_mw,
        )
        po.dblcol(self.outfile, "p_fw_alpha_mw", 0.0e0, physics_variables.p_fw_alpha_mw)
        po.dblcol(
            self.outfile, "p_fw_rad_total_mw", 0.0e0, fwbs_variables.p_fw_rad_total_mw
        )
        po.dblcol(
            self.outfile,
            "p_fw_coolant_pump_mw",
            0.0e0,
            heat_transport_variables.p_fw_coolant_pump_mw,
        )

        primsum = (
            primsum
            + fwbs_variables.p_fw_nuclear_heat_total_mw
            + physics_variables.p_fw_alpha_mw
            + fwbs_variables.p_fw_rad_total_mw
            + heat_transport_variables.p_fw_coolant_pump_mw
        )
        secsum = secsum

        po.oblnkl(self.outfile)

        po.write(self.outfile, "Blanket:")
        po.dblcol(
            self.outfile,
            "p_blkt_nuclear_heat_total_mw",
            0.0e0,
            fwbs_variables.p_blkt_nuclear_heat_total_mw,
        )
        po.write(self.outfile, "0.0e0 0.0e0 0.0e0")
        po.write(self.outfile, "0.0e0 0.0e0 0.0e0")
        po.dblcol(
            self.outfile,
            "p_blkt_coolant_pump_mw",
            0.0e0,
            heat_transport_variables.p_blkt_coolant_pump_mw,
        )

        primsum = (
            primsum
            + fwbs_variables.p_blkt_nuclear_heat_total_mw
            + heat_transport_variables.p_blkt_coolant_pump_mw
        )
        secsum = secsum

        po.oblnkl(self.outfile)

        po.write(self.outfile, "Shield:")
        po.write(
            self.outfile,
            (
                f"{fwbs_variables.p_shld_nuclear_heat_mw * heat_transport_variables.iprimshld} {fwbs_variables.p_shld_nuclear_heat_mw * (1 - heat_transport_variables.iprimshld)} {fwbs_variables.p_shld_nuclear_heat_mw}"
            ),
        )
        po.write(self.outfile, "0.0e0 0.0e0 0.0e0")
        po.write(self.outfile, "0.0e0 0.0e0 0.0e0")
        po.write(
            self.outfile,
            (
                f"{heat_transport_variables.p_shld_coolant_pump_mw * heat_transport_variables.iprimshld} {heat_transport_variables.p_shld_coolant_pump_mw * (1 - heat_transport_variables.iprimshld)} {heat_transport_variables.p_shld_coolant_pump_mw}"
            ),
        )

        primsum = (
            primsum
            + fwbs_variables.p_shld_nuclear_heat_mw * heat_transport_variables.iprimshld
            + heat_transport_variables.p_shld_coolant_pump_mw
            * heat_transport_variables.iprimshld
        )
        secsum = (
            secsum
            + fwbs_variables.p_shld_nuclear_heat_mw
            * (1 - heat_transport_variables.iprimshld)
            + heat_transport_variables.p_shld_coolant_pump_mw
            * (1 - heat_transport_variables.iprimshld)
        )

        po.oblnkl(self.outfile)

        po.write(self.outfile, "Divertor:")
        po.write(
            self.outfile,
            (
                f"{fwbs_variables.p_div_nuclear_heat_total_mw * self.iprimdiv} {fwbs_variables.p_div_nuclear_heat_total_mw * (1 - self.iprimdiv)} {fwbs_variables.p_div_nuclear_heat_total_mw}"
            ),
        )
        po.write(
            self.outfile,
            (
                f"{physics_variables.p_plasma_separatrix_mw * self.iprimdiv} {physics_variables.p_plasma_separatrix_mw * (1 - self.iprimdiv)} {physics_variables.p_plasma_separatrix_mw}"
            ),
        )
        po.write(
            self.outfile,
            (
                f"{fwbs_variables.p_div_rad_total_mw * self.iprimdiv} {fwbs_variables.p_div_rad_total_mw * (1 - self.iprimdiv)} {fwbs_variables.p_div_rad_total_mw}"
            ),
        )
        po.write(
            self.outfile,
            (
                f"{heat_transport_variables.p_div_coolant_pump_mw * self.iprimdiv} {heat_transport_variables.p_div_coolant_pump_mw * (1 - self.iprimdiv)} {heat_transport_variables.p_div_coolant_pump_mw}"
            ),
        )

        primsum = (
            primsum
            + fwbs_variables.p_div_nuclear_heat_total_mw * self.iprimdiv
            + physics_variables.p_plasma_separatrix_mw * self.iprimdiv
            + fwbs_variables.p_div_rad_total_mw * self.iprimdiv
            + heat_transport_variables.p_div_coolant_pump_mw * self.iprimdiv
        )
        secsum = (
            secsum
            + fwbs_variables.p_div_nuclear_heat_total_mw * (1 - self.iprimdiv)
            + physics_variables.p_plasma_separatrix_mw * (1 - self.iprimdiv)
            + fwbs_variables.p_div_rad_total_mw * (1 - self.iprimdiv)
            + heat_transport_variables.p_div_coolant_pump_mw * (1 - self.iprimdiv)
        )

        if physics_variables.itart == 1:
            po.oblnkl(self.outfile)
            po.write(self.outfile, "TART centrepost:")
            po.dblcol(self.outfile, "pnuc_cp", 0.0e0, fwbs_variables.pnuc_cp)
            po.write(self.outfile, "0.0e0 0.0e0 0.0e0")
            po.write(self.outfile, "0.0e0 0.0e0 0.0e0")
            po.dblcol(
                self.outfile,
                "p_cp_coolant_pump_elec_mw",
                0.0e0,
                self.p_cp_coolant_pump_elec_mw,
            )  # check

        primsum = primsum
        secsum = secsum + fwbs_variables.pnuc_cp + self.p_cp_coolant_pump_elec_mw

        po.oblnkl(self.outfile)
        po.write(self.outfile, "TF coil:")
        po.dblcol(self.outfile, "ptfnuc", 0.0e0, fwbs_variables.ptfnuc)
        po.write(self.outfile, "0.0e0 0.0e0 0.0e0")
        po.write(self.outfile, "0.0e0 0.0e0 0.0e0")
        po.write(self.outfile, "0.0e0 0.0e0 0.0e0")

        primsum = primsum
        secsum = secsum + fwbs_variables.ptfnuc

        po.oblnkl(self.outfile)
        po.write(self.outfile, "Losses to H/CD apparatus + diagnostics:")
        po.dblcol(
            self.outfile,
            "p_fw_hcd_nuclear_heat_mw",
            0.0e0,
            fwbs_variables.p_fw_hcd_nuclear_heat_mw,
        )
        po.write(self.outfile, "0.0e0 0.0e0 0.0e0")
        po.dblcol(
            self.outfile,
            "p_fw_hcd_rad_total_mw",
            0.0e0,
            fwbs_variables.p_fw_hcd_rad_total_mw,
        )
        po.write(self.outfile, "0.0e0 0.0e0 0.0e0")

        primsum = primsum
        secsum = (
            secsum
            + fwbs_variables.p_fw_hcd_nuclear_heat_mw
            + fwbs_variables.p_fw_hcd_rad_total_mw
        )

        po.oblnkl(self.outfile)
        #     write(self.outfile,'(t10,a)') repeat('-',88)
        po.write(self.outfile, (f"{primsum} {secsum} {primsum + secsum}"))
        # 10    format(t32,'neutrons',t50,f8.2,t70,f8.2,t90,f8.2)
        # 20    format(t14,'charged particle transport',t50,f8.2,t70,f8.2,t90,f8.2)
        # 30    format(t31,'radiation',t50,f8.2,t70,f8.2,t90,f8.2)
        # 40    format(t25,'coolant pumping',t50,f8.2,t70,f8.2,t90,f8.2)
        # 50    format(t34,'Totals',t50,f8.2,t70,f8.2,t90,f8.2)

        po.oblnkl(self.outfile)
        po.ovarrf(
            self.outfile,
            "Total power leaving reactor (across vacuum vessel boundary) (MW)",
            "",
            primsum + secsum + fwbs_variables.ptfnuc,
            "OP ",
        )

        po.osubhd(self.outfile, "Other secondary thermal power constituents :")
        po.ovarrf(
            self.outfile,
            "Heat removal from cryogenic plant (MW)",
            "(p_cryo_plant_electric_mw)",
            heat_transport_variables.p_cryo_plant_electric_mw,
            "OP ",
        )
        po.ovarrf(
            self.outfile,
            "Heat removal from facilities (MW)",
            "(fachtmw)",
            heat_transport_variables.fachtmw,
            "OP ",
        )

        po.ovarrf(
            self.outfile,
            "Coolant pumping efficiency losses (MW)",
            "(p_coolant_pump_loss_total_mw)",
            heat_transport_variables.p_coolant_pump_loss_total_mw,
            "OP ",
        )
        po.ovarrf(
            self.outfile,
            "Heat removal from injection power (MW)",
            "(pinjht)",
            heat_transport_variables.pinjht,
            "OP ",
        )
        po.ovarrf(
            self.outfile,
            "Heat removal from tritium plant (MW)",
            "(p_tritium_plant_electric_mw)",
            heat_transport_variables.p_tritium_plant_electric_mw,
            "OP ",
        )
        po.ovarrf(
            self.outfile,
            "Heat removal from vacuum pumps (MW)",
            "(vachtmw)",
            heat_transport_variables.vachtmw,
            "OP ",
        )
        po.ovarrf(
            self.outfile,
            "TF coil resistive power (MW)",
            "(tfcmw)",
            tfcoil_variables.tfcmw,
            "OP ",
        )

        po.oblnkl(self.outfile)
        po.ovarrf(
            self.outfile,
            "Total low-grade thermal power (MW)",
            "(p_plant_secondary_heat_mw)",
            heat_transport_variables.p_plant_secondary_heat_mw,
            "OP ",
        )
        po.ovarrf(
            self.outfile,
            "Total High-grade thermal power (MW)",
            "(pthermmw)",
            heat_transport_variables.pthermmw,
            "OP ",
        )

        po.oblnkl(self.outfile)
        po.ovarin(
            self.outfile,
            "Number of primary heat exchangers",
            "(nphx)",
            heat_transport_variables.nphx,
            "OP ",
        )

        if cost_variables.ireactor != 1:
            return

        # MDK start
        po.oblnkl(self.outfile)
        po.oblnkl(self.outfile)
        po.ocmmnt(self.outfile, "Power Balance across separatrix :")
        po.ocmmnt(self.outfile, "-------------------------------")
        po.ocmmnt(self.outfile, "Only energy deposited in the plasma is included here.")

        if physics_variables.i_rad_loss == 0:
            po.ocmmnt(
                self.outfile,
                "Total power loss is scaling power plus radiation (physics_variables.i_rad_loss = 0)",
            )
            po.ovarrf(
                self.outfile,
                "Transport power from scaling law (MW)",
                "(pscalingmw)",
                physics_variables.pscalingmw,
                "OP ",
            )
            po.ovarrf(
                self.outfile,
                "Total net radiation power (MW)",
                "(p_plasma_rad_mw)",
                physics_variables.p_plasma_rad_mw,
                "OP ",
            )
            total = physics_variables.pscalingmw + physics_variables.p_plasma_rad_mw
            po.ovarrf(self.outfile, "Total (MW)", "", total, "OP ")
        elif physics_variables.i_rad_loss == 1:
            po.ocmmnt(
                self.outfile,
                "Total power loss is scaling power plus core radiation only (physics_variables.i_rad_loss = 1)",
            )
            po.ovarrf(
                self.outfile,
                "Transport power from scaling law (MW)",
                "(pscalingmw)",
                physics_variables.pscalingmw,
                "OP ",
            )
            po.ovarrf(
                self.outfile,
                'Radiation power from inside "radius_plasma_core_norm" (MW)',
                "(pcoreradmw.)",
                physics_variables.p_plasma_inner_rad_mw,
                "OP ",
            )
            po.ovarrf(
                self.outfile,
                "Total (MW)",
                "",
                physics_variables.pscalingmw + physics_variables.p_plasma_inner_rad_mw,
                "OP ",
            )
            total = (
                physics_variables.pscalingmw + physics_variables.p_plasma_inner_rad_mw
            )
        elif physics_variables.i_rad_loss == 2:
            po.ocmmnt(
                self.outfile,
                "Total power loss is scaling power only (physics_variables.i_rad_loss = 2).",
            )
            po.ocmmnt(self.outfile, "This is not recommended for power plant models.")
            po.ovarrf(
                self.outfile,
                "Transport power from scaling law (MW)",
                "(pscalingmw)",
                physics_variables.pscalingmw,
                "OP ",
            )
            po.ovarrf(
                self.outfile, "Total (MW)", "", physics_variables.pscalingmw, "OP "
            )
            total = physics_variables.pscalingmw
        else:
            logger.error(
                f"{'The value of physics_variables.i_rad_loss appears to be invalid.'}"
            )
            po.ocmmnt(
                self.outfile,
                "ERROR: The value of physics_variables.i_rad_loss appears to be invalid.",
            )

        po.oblnkl(self.outfile)
        po.ovarrf(
            self.outfile,
            "Alpha power deposited in plasma (MW)",
            "(f_alpha_plasma*p_alpha_total_mw)",
            physics_variables.f_alpha_plasma * physics_variables.p_alpha_total_mw,
            "OP ",
        )
        po.ovarrf(
            self.outfile,
            "Power from charged products of DD and/or D-He3 fusion (MW)",
            "(p_non_alpha_charged_mw.)",
            physics_variables.p_non_alpha_charged_mw,
            "OP ",
        )
        po.ovarrf(
            self.outfile,
            "Ohmic heating (MW)",
            "(p_plasma_ohmic_mw.)",
            physics_variables.p_plasma_ohmic_mw,
            "OP ",
        )
        # if (physics_variables.i_plasma_ignited == 1) :
        #    po.ovarrf(self.outfile,'Total (MW)','',f_alpha_plasma*physics_variables.p_alpha_total_mw+physics_variables.p_non_alpha_charged_mw+p_plasma_ohmic_mw, 'OP ')
        #    po.oblnkl(self.outfile)
        #    if (abs(sum - (physics_variables.f_alpha_plasma*physics_variables.p_alpha_total_mw+physics_variables.p_non_alpha_charged_mw+physics_variables.p_plasma_ohmic_mw)) > 5.0e0) :
        #        write(*,*) 'WARNING: Power balance across separatrix is in error by more than 5 MW.'
        #    po.ocmmnt(self.outfile,'WARNING: Power balance across separatrix is in error by more than 5 MW.')
        #
        # else:
        po.ovarrf(
            self.outfile,
            "Injected power deposited in plasma (MW)",
            "(p_hcd_injected_total_mw)",
            pinj,
            "OP ",
        )
        po.ovarrf(
            self.outfile,
            "Total (MW)",
            "",
            physics_variables.f_alpha_plasma * physics_variables.p_alpha_total_mw
            + physics_variables.p_non_alpha_charged_mw
            + physics_variables.p_plasma_ohmic_mw
            + pinj,
            "OP ",
        )
        po.oblnkl(self.outfile)
        if (
            abs(
                total
                - (
                    physics_variables.f_alpha_plasma
                    * physics_variables.p_alpha_total_mw
                    + physics_variables.p_non_alpha_charged_mw
                    + physics_variables.p_plasma_ohmic_mw
                    + pinj
                )
            )
            > 5.0e0
        ):
            logger.warning(
                f"{'WARNING: Power balance across separatrix is in error by more than 5 MW.'}"
            )
            po.ocmmnt(
                self.outfile,
                "WARNING: Power balance across separatrix is in error by more than 5 MW.",
            )

        #

        po.ocmmnt(self.outfile, "Power Balance for Reactor - Summary :")
        po.ocmmnt(self.outfile, "-------------------------------------")
        po.ovarrf(
            self.outfile,
            "Fusion power (MW)",
            "(p_fusion_total_mw)",
            physics_variables.p_fusion_total_mw,
            "OP ",
        )
        po.ovarrf(
            self.outfile,
            "Power from energy multiplication in blanket and shield (MW)",
            "(emultmw)",
            fwbs_variables.emultmw,
            "OP ",
        )
        po.ovarrf(
            self.outfile,
            "Injected power (MW)",
            "(p_hcd_injected_total_mw.)",
            pinj,
            "OP ",
        )
        po.ovarrf(
            self.outfile,
            "Ohmic power (MW)",
            "(p_plasma_ohmic_mw.)",
            physics_variables.p_plasma_ohmic_mw,
            "OP ",
        )
        po.ovarrf(
            self.outfile,
            "Power deposited in primary coolant by pump (MW)",
            "(p_coolant_pump_total_mw)",
            self.p_coolant_pump_total_mw,
            "OP ",
        )
        total = (
            physics_variables.p_fusion_total_mw
            + fwbs_variables.emultmw
            + pinj
            + self.p_coolant_pump_total_mw
            + physics_variables.p_plasma_ohmic_mw
        )
        po.ovarrf(self.outfile, "Total (MW)", "", total, "OP ")
        po.oblnkl(self.outfile)
        # po.ovarrf(self.outfile,'Heat extracted from armour and first wall (MW)','(p_fw_heat_deposited_mw)',p_fw_heat_deposited_mw, 'OP ')
        po.ovarrf(
            self.outfile,
            "Heat extracted from first wall and blanket (MW)",
            "(p_fw_blkt_heat_deposited_mw)",
            self.p_fw_blkt_heat_deposited_mw,
            "OP ",
        )
        po.ovarrf(
            self.outfile,
            "Heat extracted from shield  (MW)",
            "(p_shld_heat_deposited_mw)",
            self.p_shld_heat_deposited_mw,
            "OP ",
        )
        po.ovarrf(
            self.outfile,
            "Heat extracted from divertor (MW)",
            "(p_div_heat_deposited_mw)",
            self.p_div_heat_deposited_mw,
            "OP ",
        )
        po.ovarrf(
            self.outfile,
            "Nuclear and photon power lost to H/CD system (MW)",
            "(psechcd)",
            heat_transport_variables.psechcd,
            "OP ",
        )
        po.ovarrf(
            self.outfile,
            "Nuclear power lost to TF (MW)",
            "(ptfnuc)",
            fwbs_variables.ptfnuc,
            "OP ",
        )
        po.ovarrf(
            self.outfile,
            "Total (MW)",
            "",
            self.p_fw_blkt_heat_deposited_mw
            + self.p_shld_heat_deposited_mw
            + self.p_div_heat_deposited_mw
            + heat_transport_variables.psechcd
            + fwbs_variables.ptfnuc,
            "OP ",
        )
        po.oblnkl(self.outfile)
        if (
            abs(
                total
                - (
                    self.p_fw_blkt_heat_deposited_mw
                    + self.p_shld_heat_deposited_mw
                    + self.p_div_heat_deposited_mw
                    + heat_transport_variables.psechcd
                    + fwbs_variables.ptfnuc
                )
            )
            > 5.0e0
        ):
            logger.warning(
                f"{'WARNING: Power balance for reactor is in error by more than 5 MW.'}"
            )
            po.ocmmnt(
                self.outfile,
                "WARNING: Power balance for reactor is in error by more than 5 MW.",
            )

        # Heat rejected by main power conversion circuit
        if (
            fwbs_variables.i_blkt_dual_coolant > 0
            and fwbs_variables.i_coolant_pumping == 2
        ):
            self.p_turbine_loss_mw = (
                heat_transport_variables.pthermmw - self.pthermblkt_liq
            ) * (1 - heat_transport_variables.eta_turbine) + self.pthermblkt_liq * (
                1 - heat_transport_variables.etath_liq
            )
        else:
            self.p_turbine_loss_mw = heat_transport_variables.pthermmw * (
                1 - heat_transport_variables.eta_turbine
            )

        po.ocmmnt(self.outfile, "Electrical Power Balance :")
        po.ocmmnt(self.outfile, "--------------------------")
        po.ovarrf(
            self.outfile,
            "Net electric power output(MW)",
            "(p_plant_electric_net_mw.)",
            heat_transport_variables.p_plant_electric_net_mw,
            "OP ",
        )
        po.ovarrf(
            self.outfile,
            "Required Net electric power output(MW)",
            "(pnetelin)",
            constraint_variables.pnetelin,
        )
        po.ovarrf(
            self.outfile,
            "Electric power for heating and current drive (MW)",
            "(p_hcd_electric_total_mw)",
            heat_transport_variables.p_hcd_electric_total_mw,
            "OP ",
        )
        po.ovarrf(
            self.outfile,
            "Electric power for primary coolant pumps (MW)",
            "(p_coolant_pump_elec_total_mw)",
            heat_transport_variables.p_coolant_pump_elec_total_mw,
            "OP ",
        )
        po.ovarrf(
            self.outfile,
            "Electric power for vacuum pumps (MW)",
            "(vachtmw)",
            heat_transport_variables.vachtmw,
        )
        po.ovarrf(
            self.outfile,
            "Electric power for tritium plant (MW)",
            "(p_tritium_plant_electric_mw)",
            heat_transport_variables.p_tritium_plant_electric_mw,
        )
        po.ovarrf(
            self.outfile,
            "Electric power for cryoplant (MW)",
            "(p_cryo_plant_electric_mw)",
            heat_transport_variables.p_cryo_plant_electric_mw,
            "OP ",
        )
        po.ovarrf(
            self.outfile,
            "Electric power for TF coils (MW)",
            "(tfacpd)",
            heat_transport_variables.tfacpd,
            "OP ",
        )
        po.ovarrf(
            self.outfile,
            "Electric power for PF coils (MW)",
            "(pfwpmw)",
            pfcoil_variables.pfwpmw,
            "OP ",
        )
        po.ovarrf(
            self.outfile,
            "All other internal electric power requirements (MW)",
            "(fachtmw)",
            heat_transport_variables.fachtmw,
            "OP ",
        )
        total_plant_power = (
            heat_transport_variables.p_plant_electric_net_mw
            + heat_transport_variables.p_hcd_electric_total_mw
            + heat_transport_variables.p_coolant_pump_elec_total_mw
            + heat_transport_variables.vachtmw
            + heat_transport_variables.p_tritium_plant_electric_mw
            + heat_transport_variables.p_cryo_plant_electric_mw
            + heat_transport_variables.tfacpd
            + heat_transport_variables.fachtmw
            + pfcoil_variables.pfwpmw
        )
        po.ovarrf(
            self.outfile, "Total (MW)", "(tot_plant_power)", total_plant_power, "OP "
        )
        po.ovarrf(self.outfile, "Total (MW)", "", total_plant_power, "OP ")
        po.oblnkl(self.outfile)
        po.ovarrf(
            self.outfile,
            "Gross electrical output* (MW)",
            "(p_plant_electric_gross_mw)",
            heat_transport_variables.p_plant_electric_gross_mw,
            "OP ",
        )
        po.ocmmnt(
            self.outfile, "(*Power for pumps in secondary circuit already subtracted)"
        )
        po.oblnkl(self.outfile)
        if (
            abs(total_plant_power - heat_transport_variables.p_plant_electric_gross_mw)
            > 5.0e0
        ):
            logger.warning(
                f"{'WARNING: Electrical Power balance is in error by more than 5 MW.'}"
            )
            po.ocmmnt(
                self.outfile,
                "WARNING: Electrical Power balance is in error by more than 5 MW.",
            )

        po.ocmmnt(self.outfile, "Power balance for power plant :")
        po.ocmmnt(self.outfile, "-------------------------------")
        po.ovarrf(
            self.outfile,
            "Fusion power (MW)",
            "(p_fusion_total_mw)",
            physics_variables.p_fusion_total_mw,
            "OP ",
        )
        po.ovarrf(
            self.outfile,
            "Power from energy multiplication in blanket and shield (MW)",
            "(emultmw)",
            fwbs_variables.emultmw,
            "OP ",
        )
        total_power = physics_variables.p_fusion_total_mw + fwbs_variables.emultmw
        po.ovarrf(self.outfile, "Total (MW)", "", total_power, "OP ")
        po.oblnkl(self.outfile)
        po.ovarrf(
            self.outfile,
            "Net electrical output (MW)	",
            "(p_plant_electric_net_mw)",
            heat_transport_variables.p_plant_electric_net_mw,
            "OP ",
        )
        po.ovarrf(
            self.outfile,
            "Heat rejected by main power conversion circuit (MW)",
            "(p_turbine_loss_mw)",
            self.p_turbine_loss_mw,
            "OP ",
        )
        po.ovarrf(
            self.outfile,
            "Heat rejected by other cooling circuits (MW)",
            "(p_plant_secondary_heat_mw)",
            heat_transport_variables.p_plant_secondary_heat_mw,
            "OP ",
        )
        po.ovarrf(
            self.outfile,
            "Total (MW)",
            "",
            heat_transport_variables.p_plant_electric_net_mw
            + self.p_turbine_loss_mw
            + heat_transport_variables.p_plant_secondary_heat_mw,
            "OP ",
        )
        po.oblnkl(self.outfile)
        if (
            abs(
                total_power
                - (
                    heat_transport_variables.p_plant_electric_net_mw
                    + self.p_turbine_loss_mw
                    + heat_transport_variables.p_plant_secondary_heat_mw
                )
            )
            > 5.0e0
        ):
            logger.warning(
                f"{'WARNING: Power balance for power plant is in error by more than 5 MW.'}"
            )
            po.ocmmnt(
                self.outfile,
                "WARNING: Power balance for power plant is in error by more than 5 MW.",
            )

        po.osubhd(self.outfile, "Plant efficiency measures :")
        po.ovarrf(
            self.outfile,
            "Net electric power / total nuclear power (%)",
            "(p_plant_electric_net_mw/(p_fusion_total_mw+emultmw)",
            100.0e0
            * heat_transport_variables.p_plant_electric_net_mw
            / (physics_variables.p_fusion_total_mw + fwbs_variables.emultmw),
            "OP ",
        )
        po.ovarrf(
            self.outfile,
            "Net electric power / total fusion power (%)",
            "(p_plant_electric_net_mw/p_fusion_total_mw)",
            100.0e0
            * heat_transport_variables.p_plant_electric_net_mw
            / physics_variables.p_fusion_total_mw,
            "OP ",
        )
        po.ovarrf(
            self.outfile,
            "Gross electric power* / high grade heat (%)",
            "(eta_turbine)",
            100.0e0 * heat_transport_variables.eta_turbine,
        )
        po.ocmmnt(
            self.outfile, "(*Power for pumps in secondary circuit already subtracted)"
        )
        po.ovarrf(
            self.outfile, "Recirculating power fraction", "(cirpowfr)", cirpowfr, "OP "
        )

    def power3(self, output: bool):
        """
        Calculates the time-dependent power requirements
        author: J Morris, CCFE, Culham Science Centre
        outfile : input integer : output file unit
        iprint : input integer : switch for writing to output (1=yes)
        This routine calculates the time dependent power requirements
        and outputs them to the output file
        None
        """
        p_cooling = np.zeros((6,))
        p_cryo = np.zeros((6,))
        p_vac = np.zeros((6,))
        p_tritium = np.zeros((6,))
        p_fac = np.zeros((6,))
        p_tf = np.zeros((6,))
        p_hcd = np.zeros((6,))
        p_pf = np.zeros((6,))
        p_int_tot = np.zeros((6,))
        p_gross = np.zeros((6,))

        t_cs = times_variables.t_precharge

        # Plasma current ramp up time (s)
        t_ip_up = times_variables.t_current_ramp_up

        # Plasma heating phase (s)
        t_heat = times_variables.t_fusion_ramp

        # Flat-top phase (s)
        t_flat_top = times_variables.t_burn

        # Plasma current ramp down time (s)
        t_ip_down = times_variables.t_ramp_down

        # Extra time between pulses (s)
        t_extra = times_variables.t_between_pulse

        # Continuous power usage

        # Primary pumping electrical power [MWe]
        p_cooling[0:6] = heat_transport_variables.p_coolant_pump_elec_total_mw

        # Cryoplant electrical power [MWe]
        p_cryo[0:6] = heat_transport_variables.p_cryo_plant_electric_mw

        # Vacuum electrical power [MWe]
        p_vac[0:6] = heat_transport_variables.vachtmw

        # Tritium system electrical power [MWe]
        p_tritium[0:6] = heat_transport_variables.p_tritium_plant_electric_mw

        # Facilities electrical power [MWe]
        p_fac[0:6] = heat_transport_variables.fachtmw

        # TF coil electrical power [MWe]
        p_tf[0:6] = heat_transport_variables.tfacpd

        # Total continuous power [MWe]
        p_cont_tot = p_cooling + p_cryo + p_vac + p_tritium + p_fac + p_tf

        # Intermittent power usage

        # Heating and current drive electrical power [MWe]
        p_hcd[0] = 0.0e0
        p_hcd[1] = (
            heat_transport_variables.pinjmax
            / current_drive_variables.eta_hcd_primary_injector_wall_plug
        )
        p_hcd[2] = (
            heat_transport_variables.pinjmax
            / current_drive_variables.eta_hcd_primary_injector_wall_plug
        )
        p_hcd[3] = heat_transport_variables.p_hcd_electric_total_mw
        p_hcd[4] = (
            heat_transport_variables.pinjmax
            / current_drive_variables.eta_hcd_primary_injector_wall_plug
        )
        p_hcd[5] = 0.0e0

        # PF coils electrical power [MWe]
        p_pf[0] = pf_power_variables.poloidalpower[0] / 1.0e6
        p_pf[1] = pf_power_variables.poloidalpower[1] / 1.0e6
        p_pf[2] = pf_power_variables.poloidalpower[2] / 1.0e6
        p_pf[3] = pf_power_variables.poloidalpower[3] / 1.0e6
        p_pf[4] = pf_power_variables.poloidalpower[4] / 1.0e6
        p_pf[5] = 0.0e0

        # Total intermittent power [MWe]
        p_int_tot[0:6] = p_pf + p_hcd

        # Gross power [MWe]
        p_gross[0:3] = 0.0e0
        p_gross[3] = heat_transport_variables.p_plant_electric_gross_mw
        p_gross[4:6] = 0.0e0

        # Net electric power [MWe]
        p_net = p_gross - (
            p_cooling + p_cryo + p_vac + p_fac + p_tritium + p_tf + p_pf + p_hcd
        )

        # Net electric power average [MWe]
        p_net_avg = (  # noqa: F841
            (p_net[0] * t_cs)
            + (p_net[1] * t_ip_up)
            + (p_net[2] * t_heat)
            + (p_net[3] * t_flat_top)
            + (p_net[4] * t_ip_down)
            + (p_net[5] * t_extra)
        ) / (t_cs + t_ip_up + t_heat + t_flat_top + t_ip_down + t_extra)

        # Output
        if output == 0:
            return

        po.osubhd(self.outfile, "Time-dependent power usage")

        po.write(self.outfile, "Pulse timings [s]:")
        po.oblnkl(self.outfile)
        po.write(
            self.outfile,
            "t_precharge t_current_ramp_up t_fusion_ramp t_burn t_ramp_down t_between_pulse",
        )
        po.write(self.outfile, "----- ---- ----- ----- ----- ------")
        po.write(
            self.outfile,
            (f"Duration {t_cs} {t_ip_up} {t_heat} {t_flat_top} {t_ip_down} {t_extra}"),
        )
        po.write(self.outfile, "------ ----- ---- ----- ----- ----- ------")
        po.oblnkl(self.outfile)

        po.write(self.outfile, "Continous power usage [MWe]:")
        po.oblnkl(self.outfile)
        po.write(
            self.outfile,
            "System t_precharge t_current_ramp_up t_fusion_ramp t_burn t_ramp_down t_between_pulse",
        )
        po.write(self.outfile, "------ ----- ---- ----- ----- ----- ------")
        po.write(
            self.outfile,
            (
                f"Primary cooling {p_cooling[0]} {p_cooling[1]} {p_cooling[2]} {p_cooling[3]} {p_cooling[4]} {p_cooling[5]}"
            ),
        )
        po.write(
            self.outfile,
            (
                f"Cryoplant {p_cryo[0]} {p_cryo[1]} {p_cryo[2]} {p_cryo[3]} {p_cryo[4]} {p_cryo[5]}"
            ),
        )
        po.write(
            self.outfile,
            (
                f"Vacuum {p_vac[0]} {p_vac[1]} {p_vac[2]} {p_vac[3]} {p_vac[4]} {p_vac[5]}"
            ),
        )
        po.write(
            self.outfile,
            (
                f"Tritium {p_tritium[0]} {p_tritium[1]} {p_tritium[2]} {p_tritium[3]} {p_tritium[4]} {p_tritium[5]}"
            ),
        )
        po.write(
            self.outfile,
            (f"TF {p_tf[0]} {p_tf[1]} {p_tf[2]} {p_tf[3]} {p_tf[4]} {p_tf[5]}"),
        )
        po.write(
            self.outfile,
            (
                f"Facilities {p_fac[0]} {p_fac[1]} {p_fac[2]} {p_fac[3]} {p_fac[4]} {p_fac[5]}"
            ),
        )

        po.write(self.outfile, "------ ----- ---- ----- ----- ----- ------")
        po.write(
            self.outfile,
            (
                f"Total {p_cont_tot[0]} {p_cont_tot[1]} {p_cont_tot[2]} {p_cont_tot[3]} {p_cont_tot[4]} {p_cont_tot[5]}"
            ),
        )
        po.write(self.outfile, "------ ----- ---- ----- ----- ----- ------")
        po.oblnkl(self.outfile)

        po.write(self.outfile, "Intermittent power usage [MWe]:")
        po.oblnkl(self.outfile)
        po.write(
            self.outfile,
            "System t_precharge t_current_ramp_up t_fusion_ramp t_burn t_ramp_down t_between_pulse",
        )
        po.write(self.outfile, "------ ----- ---- ----- ----- ----- ------")
        po.write(
            self.outfile,
            (
                f"H & CD {p_hcd[0]} {p_hcd[1]} {p_hcd[2]} {p_hcd[3]} {p_hcd[4]} {p_hcd[5]}"
            ),
        )
        po.write(
            self.outfile,
            (f"PF {p_pf[0]} {p_pf[1]} {p_pf[2]} {p_pf[3]} {p_pf[4]} {p_pf[5]}"),
        )

        po.write(self.outfile, "------ ----- ---- ----- ----- ----- ------")
        po.write(
            self.outfile,
            (
                f"Total {p_int_tot[0]} {p_int_tot[1]} {p_int_tot[2]} {p_int_tot[3]} {p_int_tot[4]} {p_int_tot[5]}"
            ),
        )
        po.write(self.outfile, "------ ----- ---- ----- ----- ----- ------")

        po.oblnkl(self.outfile)

        po.write(self.outfile, "Power production [MWe]:")
        po.oblnkl(self.outfile)
        po.write(
            self.outfile,
            " t_precharge t_current_ramp_up t_fusion_ramp t_burn t_ramp_down t_between_pulse avg",
        )
        po.write(self.outfile, " ----- ---- ----- ----- ----- ------ ---")
        po.write(
            self.outfile,
            (
                f"Gross power {p_gross[0]} {p_gross[1]} {p_gross[2]} {p_gross[3]} {p_gross[4]} {p_gross[5]}"
            ),
        )
        po.write(
            self.outfile,
            (
                f"Net power {p_net[0]} {p_net[1]} {p_net[2]} {p_net[3]} {p_net[4]} {p_net[5]}"
            ),
        )
        po.write(self.outfile, "------ ----- ---- ----- ----- ----- ------")

        po.oblnkl(self.outfile)

    # 10    format(t20,a20,t40,a8,t50,a8,t60,a8,t70,a8,t80,a8,t90,a8)
    # 20    format(t20,a20,t40,f8.2,t50,f8.2,t60,f8.2,t70,f8.2,t80,f8.2,t90,f8.2,t100,f8.2)
    # 30    format(t20,a20,t40,a8,t50,a8,t60,a8,t70,a8,t80,a8,t90,a8,t100,a8)
    # 40    format(t20,a20,t40,f8.2,t50,f8.2,t60,f8.2,t70,f8.2,t80,f8.2,t90,f8.2,t100,f8.2,t110,f8.2)

    def cryo(
        self,
        i_tf_sup,
        tfcryoarea,
        coldmass,
        ptfnuc,
        ensxpfm,
        t_pulse_repetition,
        cpttf,
        n_tf_coils,
    ):
        """
        Calculates cryogenic loads
        author: P J Knight, CCFE, Culham Science Centre
        itfsup : input integer : Switch denoting whether TF coils are
        superconducting
        tfcryoarea : input real : Surface area of toroidal shells covering TF coils (m2)
        coldmass : input real : Mass of cold (cryogenic) components (kg),
        including TF coils, PF coils, cryostat, and
        intercoil structure
        ptfnuc : input real : Nuclear heating in TF coils (MW)
        ensxpfm : input real : Maximum PF coil stored energy (MJ)
        t_pulse_repetition : input real : Pulse length of cycle (s)
        cpttf : input real : Current per turn in TF coils (A)
        tfno : input real : Number of TF coils
        helpow : output real : Helium heat removal at cryo temperatures (W)
        This routine calculates the cryogenic heat load.
        D. Slack memo SCMDG 88-5-1-059, LLNL ITER-88-054, Aug. 1988
        """
        self.qss = 4.3e-4 * coldmass
        if i_tf_sup == 1:
            self.qss = self.qss + 2.0e0 * tfcryoarea

        #  Nuclear heating of TF coils (W) (zero if resistive)
        if fwbs_variables.inuclear == 0 and i_tf_sup == 1:
            fwbs_variables.qnuc = 1.0e6 * ptfnuc
        # Issue #511: if fwbs_variables.inuclear = 1 : fwbs_variables.qnuc is input.

        #  AC losses
        self.qac = 1.0e3 * ensxpfm / t_pulse_repetition

        #  Current leads
        if i_tf_sup == 1:
            self.qcl = 13.6e-3 * n_tf_coils * cpttf
        else:
            self.qcl = 0.0e0

        #  45% extra miscellaneous, piping and reserves
        self.qmisc = 0.45e0 * (self.qss + fwbs_variables.qnuc + self.qac + self.qcl)
        return max(
            0.0e0,
            self.qmisc + self.qss + fwbs_variables.qnuc + self.qac + self.qcl,
        )

    def plant_thermal_efficiency(self, eta_turbine):
        """
        Calculates the thermal efficiency of the power conversion cycle
        author: P J Knight, CCFE, Culham Science Centre
        author: C Harrington, CCFE, Culham Science Centre
        eta_turbine : input/output real : thermal to electric conversion efficiency
        This routine calculates the thermal efficiency of the power conversion cycle.
        This gives the gross power of the plant, i.e. the primary coolant pumping
        power is not subtracted at this point; however, the pumping of the
        secondary coolant is accounted for.
        <P>If i_thermal_electric_conversion = 0, 1, a set efficiency for the chosen blanket design is used,
        taken from cycle modelling studies.
        <P>If i_thermal_electric_conversion > 1, the outlet temperature from the first wall
        and breeder zone is used to calculate an efficiency, using a simple relationship
        between eta_turbine and temp_blkt_coolant_out again obtained from previous studies.
        C. Harrington, K:Power Plant Physics and Technology  PROCESS  blanket_model
         New Power Module Harrington  Cycle correlations  Cycle correlations.xls
        """
        if fwbs_variables.i_thermal_electric_conversion == 0:
            #  CCFE HCPB Model (with or without TBR)
            if fwbs_variables.i_blanket_type == 1:
                #  HCPB, efficiency taken from M. Kovari 2016
                # "PROCESS": A systems code for fusion power plants - Part 2: Engineering
                # https://www.sciencedirect.com/science/article/pii/S0920379616300072
                # Feedheat & reheat cycle assumed
                eta_turbine = 0.411e0
            else:
                logger.log(f"{'i_blanket_type does not have a value in range 1-3.'}")

            #  Etath from reference. Div power to primary
        elif fwbs_variables.i_thermal_electric_conversion == 1:
            #  CCFE HCPB Model (with or without TBR)
            if fwbs_variables.i_blanket_type == 1:
                #  HCPB, efficiency taken from M. Kovari 2016
                # "PROCESS": A systems code for fusion power plants - Part 2: Engineering
                # https://www.sciencedirect.com/science/article/pii/S0920379616300072
                # Feedheat & reheat cycle assumed
                eta_turbine = 0.411e0 - self.delta_eta
            else:
                logger.log(f"{'i_blanket_type does not have a value in range 1-3.'}")

            #  User input used, eta_turbine not changed
        elif fwbs_variables.i_thermal_electric_conversion == 2:
            return eta_turbine
            # Do nothing

            #  Steam Rankine cycle to be used
        elif fwbs_variables.i_thermal_electric_conversion == 3:
            #  CCFE HCPB Model (with or without TBR)
            if fwbs_variables.i_blanket_type == 1:
                #  If coolant is helium, the steam cycle is assumed to be superheated
                #  and a different correlation is used. The turbine inlet temperature
                #  is assumed to be 20 degrees below the primary coolant outlet
                #  temperature, as was stated for steam rankine cycle for Helium in
                #  M. Kovari 2016, "PROCESS": A systems code for fusion power plants
                #  - Part 2: Engineering
                #  https://www.sciencedirect.com/science/article/pii/S0920379616300072

                #  Superheated steam Rankine cycle correlation (C. Harrington)
                #  Range of validity: 657 K < heat_transport_variables.temp_turbine_coolant_in < 915 K
                heat_transport_variables.temp_turbine_coolant_in = (
                    fwbs_variables.temp_blkt_coolant_out - 20.0e0
                )
                if (heat_transport_variables.temp_turbine_coolant_in < 657.0e0) or (
                    heat_transport_variables.temp_turbine_coolant_in > 915.0e0
                ):
                    error_handling.idiags[0] = 2
                    error_handling.fdiags[0] = (
                        heat_transport_variables.temp_turbine_coolant_in
                    )
                    error_handling.report_error(166)

                eta_turbine = (
                    0.1802e0 * np.log(heat_transport_variables.temp_turbine_coolant_in)
                    - 0.7823
                    - self.delta_eta
                )

            else:
                logger.log(f"{'i_blanket_type does not have a value in range 1-3.'}")

            #  Supercritical CO2 cycle to be used
        elif fwbs_variables.i_thermal_electric_conversion == 4:
            #  The same temperature/efficiency correlation is used regardless of
            #  primary coolant choice.  The turbine inlet temperature is assumed to
            #  be 20 degrees below the primary coolant outlet temperature.
            #  s-CO2 can in theory be used for both helium and water primary coolants
            #  so no differentiation is made, but for water the efficiency will be
            #  very low and the correlation will reflect this.

            #  Supercritical CO2 cycle correlation (C. Harrington)
            #  Range of validity: 408 K < heat_transport_variables.temp_turbine_coolant_in < 1023 K
            heat_transport_variables.temp_turbine_coolant_in = (
                fwbs_variables.temp_blkt_coolant_out - 20.0e0
            )
            if (heat_transport_variables.temp_turbine_coolant_in < 408.0e0) or (
                heat_transport_variables.temp_turbine_coolant_in > 1023.0e0
            ):
                error_handling.idiags[0] = 3
                error_handling.fdiags[0] = (
                    heat_transport_variables.temp_turbine_coolant_in
                )
                error_handling.report_error(166)

            eta_turbine = (
                0.4347e0 * np.log(heat_transport_variables.temp_turbine_coolant_in)
                - 2.5043e0
            )

        else:
            logger.log(
                f"{'i_thermal_electric_conversion does not appear to have a value within its range (0-4)'}"
            )
        return eta_turbine

    def plant_thermal_efficiency_2(self, etath_liq):
        """
        Calculates the thermal efficiency of the power conversion cycle
        for the liquid metal breeder
        """
        if fwbs_variables.secondary_cycle_liq == 2:
            #  User input used, eta_turbine not changed
            return etath_liq

        if fwbs_variables.secondary_cycle_liq == 4:
            #  Supercritical CO2 cycle to be used
            #  Supercritical CO2 cycle correlation (C. Harrington)
            #  Range of validity: 408 K < heat_transport_variables.temp_turbine_coolant_in < 1023 K
            heat_transport_variables.temp_turbine_coolant_in = (
                fwbs_variables.outlet_temp_liq - 20.0e0
            )
            if (heat_transport_variables.temp_turbine_coolant_in < 408.0e0) or (
                heat_transport_variables.temp_turbine_coolant_in > 1023.0e0
            ):
                error_handling.idiags[0] = 3
                error_handling.fdiags[0] = (
                    heat_transport_variables.temp_turbine_coolant_in
                )
                error_handling.report_error(166)

            return (
                0.4347e0 * np.log(heat_transport_variables.temp_turbine_coolant_in)
                - 2.5043e0
            )

        raise ProcessValueError(
            f"secondary_cycle_liq ={fwbs_variables.secondary_cycle_liq} is an invalid option."
        )

    def tfpwr(self, output: bool):
        """
        TF coil power supply requirements for resistive coils
        author: P J Knight, CCFE, Culham Science Centre
        outfile : input integer : output file unit
        iprint : input integer : switch for writing to output (1=yes)
        This routine calculates the power conversion requirements for
        resistive TF coils, or calls <CODE>tfpwcall</CODE> if the TF
        coils are superconducting.
        None
        """
        if tfcoil_variables.i_tf_sup != 1:
            # Cross-sectional area of bus
            # tfcoil_variables.cpttf  - current per TFC turn (A)
            # tfcoil_variables.j_tf_bus   - bus current density (A/m2)
            a_tf_bus = tfcoil_variables.cpttf / tfcoil_variables.j_tf_bus

            # Bus resistance [ohm]
            # Bus resistivity (tfcoil_variables.rho_tf_bus)
            # Issue #1253: there was a fudge here to set the bus bar resistivity equal
            # to the TF conductor resistivity. I have removed this.
            tfbusres = (
                tfcoil_variables.rho_tf_bus * tfcoil_variables.len_tf_bus / a_tf_bus
            )

            #  Bus mass (kg)
            tfcoil_variables.m_tf_bus = (
                tfcoil_variables.len_tf_bus * a_tf_bus * constants.dcopper
            )

            #  Total maximum impedance MDK actually just fixed resistance
            res_tf_system_total = (
                tfcoil_variables.n_tf_coils * tfcoil_variables.res_tf_leg
                + (tfcoil_variables.p_cp_resistive / tfcoil_variables.c_tf_total**2)
                + tfbusres
            )

            #  No reactive portion of the voltage is included here - assume long ramp times
            #  MDK This is steady state voltage, not "peak" voltage
            tfcoil_variables.vtfkv = (
                1.0e-3
                * res_tf_system_total
                * tfcoil_variables.cpttf
                / tfcoil_variables.n_tf_coils
            )

            # Resistive powers (MW):
            tfcoil_variables.tfcpmw = (
                1.0e-6 * tfcoil_variables.p_cp_resistive
            )  # inboard legs (called centrepost, CP for tart design)
            tfcoil_variables.tflegmw = (
                1.0e-6 * tfcoil_variables.p_tf_leg_resistive
            )  # outboard legs
            tfcoil_variables.tfjtsmw = 1.0e-6 * tfcoil_variables.pres_joints  # Joints
            tfbusmw = (
                1.0e-6 * tfcoil_variables.cpttf**2 * tfbusres
            )  # TF coil bus => Dodgy #

            #  TF coil reactive power
            #  Set reactive power to 0, since ramp up can be long
            #  The TF coil can be ramped up as slowly as you like
            #  (although this will affect the time to recover from a magnet quench).
            #     tfreacmw = 1.0e-6 * 1.0e9 * estotf/(t_current_ramp_up + t_precharge)
            #                                 estotf(=estotftgj/tfcoil_variables.n_tf_coils) has been removed (#199 #847)
            tfreacmw = 0.0e0

            # Total power consumption (MW)
            tfcoil_variables.tfcmw = (
                tfcoil_variables.tfcpmw
                + tfcoil_variables.tflegmw
                + tfbusmw
                + tfreacmw
                + tfcoil_variables.tfjtsmw
            )

            # Total steady state AC power demand (MW)
            heat_transport_variables.tfacpd = (
                tfcoil_variables.tfcmw / heat_transport_variables.etatf
            )

        else:  # Superconducting TF coil option
            self.tfpwcall(output)
            return

        # Output section
        if output == 0:
            return
        # Clarify that these outputs are for resistive coils only
        po.oheadr(self.outfile, "Resistive TF Coil Power Conversion")
        po.ovarre(self.outfile, "Bus resistance (ohm)", "(tfbusres)", tfbusres, "OP ")
        po.ovarre(
            self.outfile,
            "Bus current density (A/m2)",
            "(j_tf_bus)",
            tfcoil_variables.j_tf_bus,
        )
        po.ovarre(
            self.outfile,
            "Bus length - all coils (m)",
            "(len_tf_bus)",
            tfcoil_variables.len_tf_bus,
        )
        po.ovarre(
            self.outfile,
            "Bus mass (kg)",
            "(m_tf_bus)",
            tfcoil_variables.m_tf_bus,
            "OP ",
        )
        # po.ovarre(outfile,'Maximum impedance (ohm)','(ztot)',ztot)
        po.ovarre(
            self.outfile,
            "Total resistance for TF coil set (ohm)",
            "(res_tf_system_total)",
            res_tf_system_total,
            "OP ",
        )
        # po.ovarre(outfile,'Peak voltage per coil (kV)','(vtfkv)',vtfkv)
        po.ovarre(
            self.outfile,
            "Steady-state voltage per coil (kV)",
            "(vtfkv)",
            tfcoil_variables.vtfkv,
            "OP ",
        )
        # po.ovarre(outfile,'Peak power (MW)','(tfcmw..)',tfcmw)
        po.ovarre(
            self.outfile,
            "Total power dissipation in TF coil set (MW)",
            "(tfcmw..)",
            tfcoil_variables.tfcmw,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Power dissipation in TF coil set: inboard legs (MW)",
            "(tfcpmw)",
            tfcoil_variables.tfcpmw,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Power dissipation in TF coil set: outboard legs (MW)",
            "(tflegmw)",
            tfcoil_variables.tflegmw,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Power dissipation in TF coil set: buses",
            "(tfbusmw)",
            tfbusmw,
            "OP ",
        )
        if tfcoil_variables.i_cp_joints != 0:
            po.ovarre(
                self.outfile,
                "Power dissipation in TF coil set: joints",
                "(tfjtsmw)",
                tfcoil_variables.tfjtsmw,
                "OP ",
            )

        # Reactive poower has been set to zero.
        # po.ovarre(outfile,'TF coil reactive power (MW)','(tfreacmw)', tfreacmw)

    def tfpwcall(self, output: bool):
        """
        Calls the TF coil power conversion routine for
        superconducting coils
        author: P J Knight, CCFE, Culham Science Centre
        author: P C Shipe, ORNL
        outfile : input integer : output file unit
        iprint : input integer : switch for writing to output (1=yes)
        This routine calls routine <CODE>tfcpwr</CODE> to calculate
        the power conversion requirements for superconducting TF coils.
        None
        """
        ettfmj = tfcoil_variables.estotftgj / tfcoil_variables.n_tf_coils * 1.0e3

        #  TF coil current (kA)

        itfka = 1.0e-3 * tfcoil_variables.cpttf

        (
            tfcoil_variables.tfckw,
            tfcoil_variables.len_tf_bus,
            tfcoil_variables.drarea,
            buildings_variables.tfcbv,
            heat_transport_variables.tfacpd,
        ) = self.tfcpwr(
            output,
            itfka,
            physics_variables.rmajor,
            tfcoil_variables.n_tf_coils,
            tfcoil_variables.vtfskv,
            ettfmj,
            tfcoil_variables.res_tf_leg,
        )

    def tfcpwr(self, output: bool, itfka, rmajor, ntfc, vtfskv, ettfmj, rptfc):
        """
        Calculates the TF coil power conversion system parameters
        for superconducting coils
        author: P J Knight, CCFE, Culham Science Centre
        author: P C Shipe, ORNL
        This routine calculates the TF power conversion system
        parameters:  floor space, power supplies, bussing,
        coil protection equipment, and the associated controls
        and instrumentation. It was originally written by G. Gorker,
        FEDC/ORNL, April 1987, modified by J. Galambos in 1991 to
        run in TETRA, and included in PROCESS in 1992 by P. C. Shipe.
        None
        """

        ncpbkr = 1.0e0  # number of TF coils per circuit breaker
        djmka = 0.125e0  # design current density of TF bus, kA/cm2
        rtfps = 1.05e0  # rating factor for TF coil power supplies
        fspc1 = 0.15e0  # floor space coefficient for power supplies
        fspc2 = 0.8e0  # floor space coefficient for circuit breakers
        fspc3 = 0.4e0  # floor space coefficient for load centres

        if rptfc == 0.0e0:
            tchghr = 4.0e0  # charge time of the coils, hours
            nsptfc = 1.0e0  # superconducting (1.0 = superconducting, 0.0 = resistive)
        else:
            tchghr = 0.16667e0  # charge time of the coils, hours
            nsptfc = 0.0e0  # resistive (1.0 = superconducting, 0.0 = resistive)

        #  Total steady state TF coil AC power demand (summed later)
        tfacpd = 0.0e0

        #  Stored energy of all TF coils, MJ
        ettfc = ntfc * ettfmj

        #  Inductance of all TF coils, Henries
        ltfth = 2.0e0 * ettfc / itfka**2

        #  Number of circuit breakers
        ntfbkr = ntfc / ncpbkr

        #  Inductance per TF coil, Henries
        lptfcs = ltfth / ntfc

        #  Aluminium bus section area, sq cm
        albusa = itfka / djmka

        #  Total TF system bus length, m
        len_tf_bus = (
            8.0e0 * np.pi * rmajor
            + (1.0e0 + ntfbkr) * (12.0e0 * rmajor + 80.0e0)
            + 0.2e0 * itfka * np.sqrt(ntfc * rptfc * 1000.0e0)
        )

        #  Aluminium bus weight, tonnes
        albuswt = 2.7e0 * albusa * len_tf_bus / 1.0e4

        #  Total resistance of TF bus, ohms
        # rtfbus = 2.62e-4 * len_tf_bus / albusa
        rtfbus = tfcoil_variables.rho_tf_bus * len_tf_bus / (albusa / 10000)

        #  Total voltage drop across TF bus, volts
        vtfbus = 1000.0e0 * itfka * rtfbus

        #  Total resistance of the TF coils, ohms
        rcoils = ntfc * rptfc

        #  Total impedance, ohms
        ztotal = rtfbus + rcoils + ltfth / (3600.0e0 * tchghr)

        #  Charging voltage for the TF coils, volts
        tfcv = 1000.0e0 * itfka * ztotal

        #  Number of TF power modules
        ntfpm = (itfka * (1.0e0 + nsptfc)) / 5.0e0

        #  TF coil power module voltage, volts
        tfpmv = rtfps * tfcv / (1.0e0 + nsptfc)

        #  TF coil power supply voltage, volts
        tfpsv = rtfps * tfcv

        #  Power supply current, kA
        tfpska = rtfps * itfka

        #  TF power module current, kA
        tfpmka = rtfps * itfka / (ntfpm / (1.0e0 + nsptfc))

        #  TF power module power, kW
        tfpmkw = tfpmv * tfpmka

        #  Available DC power for charging the TF coils, kW
        tfckw = tfpmkw * ntfpm

        #  Peak AC power needed to charge coils, kW
        tfackw = tfckw / 0.9e0

        #  Resistance of dump resistor, ohms
        r1dump = nsptfc * vtfskv * ncpbkr / itfka

        #  Time constant, s
        ttfsec = lptfcs * ncpbkr / (r1dump * nsptfc + rptfc * (1.0e0 - nsptfc))

        #  Number of dump resistors
        ndumpr = ntfbkr * 4.0e0

        #  Peak power to a dump resistor during quench, MW
        r1ppmw = nsptfc * r1dump * (itfka / 2.0e0) ** 2

        #  Energy to dump resistor during quench, MJ
        r1emj = nsptfc * ettfc / (ndumpr + 0.0001e0)

        #  Total TF coil peak resistive power demand, MVA
        rpower = (ntfc * rptfc + rtfbus) * itfka**2

        #  Total TF coil peak inductive power demand, MVA
        xpower = ltfth / (3600.0e0 * tchghr) * itfka**2

        #  Building space:
        #  Power modules floor space, m2
        part1 = fspc1 * ntfpm * tfpmkw**0.667e0

        #  Circuit breakers floor space, m2
        part2 = fspc2 * ntfbkr * (vtfskv * itfka) ** 0.667e0

        #  Load centres floor space, m2
        part3 = (
            fspc3 * (tfackw / (2.4e0 * nsptfc + 13.8e0 * (1.0e0 - nsptfc))) ** 0.667e0
        )

        #  Power conversion building floor area, m2
        tfcfsp = part1 + part2 + part3

        #  Dump resistor floor area, m2
        drarea = 0.5e0 * ndumpr * (1.0e0 + r1emj) ** 0.667e0

        #  Total TF coil power conversion building volume, m3
        tfcbv = 6.0e0 * tfcfsp

        #  TF coil AC inductive power demand, MW
        xpwrmw = xpower / 0.9e0

        #  Total steady state AC power demand, MW
        tfacpd = tfacpd + rpower / heat_transport_variables.etatf
        #  Total TF coil power conversion building floor area, m2

        # tftsp = tfcfsp
        #  Total TF coil power conversion building volume, m3

        # tftbv = tfcbv

        #  Output section
        if output:
            po.oheadr(self.outfile, "Superconducting TF Coil Power Conversion")
            po.ovarre(self.outfile, "TF coil current (kA)", "(itfka)", itfka, "OP ")
            po.ovarre(self.outfile, "Number of TF coils", "(ntfc)", ntfc)
            po.ovarre(
                self.outfile,
                "Voltage across a TF coil during quench (kV)",
                "(vtfskv)",
                vtfskv,
                "OP ",
            )
            po.ovarre(self.outfile, "TF coil charge time (hours)", "(tchghr)", tchghr)
            po.ovarre(
                self.outfile,
                "Total inductance of TF coils (H)",
                "(ltfth)",
                ltfth,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Total resistance of TF coils (ohm)",
                "(rcoils)",
                rcoils,
                "OP ",
            )
            # MDK Remove this as it leads to confusion between (a) total inductance/n_tf_coils, or (b)
            #     self-inductance of one single coil
            # po.ovarre(outfile,'Inductance per TF coil (H)','(lptfcs)',lptfcs, 'OP ')
            po.ovarre(self.outfile, "TF coil charging voltage (V)", "(tfcv)", tfcv)
            po.ovarre(self.outfile, "Number of DC circuit breakers", "(ntfbkr)", ntfbkr)
            po.ovarre(self.outfile, "Number of dump resistors", "(ndumpr)", ndumpr)
            po.ovarre(
                self.outfile,
                "Resistance per dump resistor (ohm)",
                "(r1dump)",
                r1dump,
                "OP ",
            )
            po.ovarre(
                self.outfile, "Dump resistor peak power (MW)", "(r1ppmw)", r1ppmw, "OP "
            )
            po.ovarre(
                self.outfile,
                "Energy supplied per dump resistor (MJ)",
                "(r1emj)",
                r1emj,
                "OP ",
            )
            po.ovarre(
                self.outfile, "TF coil L/R time constant (s)", "(ttfsec)", ttfsec, "OP "
            )

            po.ovarre(self.outfile, "Power supply voltage (V)", "(tfpsv)", tfpsv, "OP ")
            po.ovarre(
                self.outfile, "Power supply current (kA)", "(tfpska)", tfpska, "OP "
            )
            po.ovarre(
                self.outfile, "DC power supply rating (kW)", "(tfckw)", tfckw, "OP "
            )
            po.ovarre(
                self.outfile, "AC power for charging (kW)", "(tfackw)", tfackw, "OP "
            )
            po.ovarre(
                self.outfile, "TF coil resistive power (MW)", "(rpower)", rpower, "OP "
            )

            po.ovarre(
                self.outfile, "TF coil inductive power (MVA)", "(xpower)", xpower, "OP "
            )
            po.ovarre(
                self.outfile, "Aluminium bus current density (kA/cm2)", "(djmka)", djmka
            )
            po.ovarre(
                self.outfile,
                "Aluminium bus cross-sectional area (cm2)",
                "(albusa)",
                albusa,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Total length of TF coil bussing (m)",
                "(len_tf_bus)",
                len_tf_bus,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Aluminium bus weight (tonnes)",
                "(albuswt)",
                albuswt,
                "OP ",
            )

            po.ovarre(
                self.outfile,
                "Total TF coil bus resistance (ohm)",
                "(rtfbus)",
                rtfbus,
                "OP ",
            )
            po.ovarre(
                self.outfile, "TF coil bus voltage drop (V)", "(vtfbus)", vtfbus, "OP "
            )
            po.ovarre(
                self.outfile, "Dump resistor floor area (m2)", "(drarea)", drarea, "OP "
            )
            po.ovarre(
                self.outfile,
                "TF coil power conversion floor space (m2)",
                "(tfcfsp)",
                tfcfsp,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "TF coil power conv. building volume (m3)",
                "(tfcbv)",
                tfcbv,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "TF coil AC inductive power demand (MW)",
                "(xpwrmw)",
                xpwrmw,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Total steady state AC power demand (MW)",
                "(tfacpd)",
                tfacpd,
                "OP ",
            )

        return (tfckw, len_tf_bus, drarea, tfcbv, tfacpd)


def init_pf_power_variables():
    """Initialise PF coil power variables"""
    pf_power_variables.acptmax = 0.0
    pf_power_variables.ensxpfm = 0.0
    pf_power_variables.iscenr = 2
    pf_power_variables.pfckts = 0.0
    pf_power_variables.spfbusl = 0.0
    pf_power_variables.spsmva = 0.0
    pf_power_variables.srcktpm = 0.0
    pf_power_variables.vpfskv = 0.0
    pf_power_variables.peakpoloidalpower = 0.0
    pf_power_variables.maxpoloidalpower = 1000.0
    pf_power_variables.poloidalpower[:] = 0.0


def init_heat_transport_variables():
    """Initialise heat transport variables"""
    heat_transport_variables.p_plant_electric_base = 5.0e6
    heat_transport_variables.p_cryo_plant_electric_mw = 0.0
    heat_transport_variables.p_cryo_plant_electric_max_mw = 50.0
    heat_transport_variables.f_crypmw = 1.0
    heat_transport_variables.etatf = 0.9
    heat_transport_variables.eta_turbine = 0.35
    heat_transport_variables.etath_liq = 0.35
    heat_transport_variables.fachtmw = 0.0
    heat_transport_variables.fcsht = 0.0
    heat_transport_variables.fgrosbop = 0.0
    heat_transport_variables.fmgdmw = 0.0
    heat_transport_variables.fpumpblkt = 0.005
    heat_transport_variables.fpumpdiv = 0.005
    heat_transport_variables.fpumpfw = 0.005
    heat_transport_variables.fpumpshld = 0.005
    heat_transport_variables.helpow = 0.0
    heat_transport_variables.helpow_cryal = 0.0
    heat_transport_variables.p_coolant_pump_elec_total_mw = 0.0
    heat_transport_variables.p_blkt_coolant_pump_mw = 0.0
    heat_transport_variables.p_blkt_breeder_pump_mw = 0.0
    heat_transport_variables.htpmw_blkt_tot = 0.0
    heat_transport_variables.p_div_coolant_pump_mw = 0.0
    heat_transport_variables.p_fw_coolant_pump_mw = 0.0
    heat_transport_variables.p_shld_coolant_pump_mw = 0.0
    heat_transport_variables.p_coolant_pump_loss_total_mw = 0.0
    heat_transport_variables.ipowerflow = 1
    heat_transport_variables.iprimshld = 1
    heat_transport_variables.nphx = 0
    heat_transport_variables.pacpmw = 0.0
    heat_transport_variables.peakmva = 0.0
    heat_transport_variables.p_fw_div_heat_deposited_mw = 0.0
    heat_transport_variables.p_plant_electric_gross_mw = 0.0
    heat_transport_variables.pinjht = 0.0
    heat_transport_variables.pinjmax = 120.0
    heat_transport_variables.p_hcd_electric_total_mw = 0.0
    heat_transport_variables.p_hcd_secondary_electric_mw = 0.0
    heat_transport_variables.p_plant_electric_net_mw = 0.0
    heat_transport_variables.p_plant_electric_recirc_mw = 0.0
    heat_transport_variables.priheat = 0.0
    heat_transport_variables.psecdiv = 0.0
    heat_transport_variables.psechcd = 0.0
    heat_transport_variables.p_plant_secondary_heat_mw = 0.0
    heat_transport_variables.pseclossmw = 0.0
    heat_transport_variables.psecshld = 0.0
    heat_transport_variables.pthermmw = 0.0
    heat_transport_variables.pwpm2 = 150.0
    heat_transport_variables.tfacpd = 0.0
    heat_transport_variables.tlvpmw = 0.0
    heat_transport_variables.p_tritium_plant_electric_mw = 15.0
    heat_transport_variables.temp_turbine_coolant_in = 0.0
    heat_transport_variables.vachtmw = 0.5
