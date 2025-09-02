import logging
import math

import numpy as np
import scipy as sp

from process import constants
from process import process_output as po
from process.data_structure import (
    build_variables,
    buildings_variables,
    cost_variables,
    current_drive_variables,
    fwbs_variables,
    heat_transport_variables,
    numerics,
    pf_power_variables,
    pfcoil_variables,
    physics_variables,
    power_variables,
    primary_pumping_variables,
    structure_variables,
    tfcoil_variables,
    times_variables,
)
from process.exceptions import ProcessValueError

logger = logging.getLogger(__name__)


class Power:
    def __init__(self):
        self.outfile = constants.NOUT
        self.mfile = constants.MFILE

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
        powpfii = np.zeros((pfcoil_variables.NGC2,))
        cktr = np.zeros((pfcoil_variables.NGC2,))
        pfcr = np.zeros((pfcoil_variables.NGC2,))
        albusa = np.zeros((pfcoil_variables.NGC2,))
        pfbusr = np.zeros((pfcoil_variables.NGC2,))
        pfcr = np.zeros((pfcoil_variables.NGC2,))
        cktr = np.zeros((pfcoil_variables.NGC2,))
        rcktvm = np.zeros((pfcoil_variables.NGC2,))
        rcktpm = np.zeros((pfcoil_variables.NGC2,))
        vpfi = np.zeros((pfcoil_variables.NGC2,))
        psmva = np.zeros((pfcoil_variables.NGC2,))
        poloidalenergy = np.zeros((6,))
        inductxcurrent = np.zeros((6,))
        pfdissipation = np.zeros((5,))

        #  Bus length
        pfbusl = 8.0e0 * physics_variables.rmajor + 140.0e0

        #  Find power requirements for PF coils at times_variables.t_pulse_cumulative(ktim)

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
        delktim = times_variables.t_plant_pulse_plasma_current_ramp_up

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

                    #  Voltage in circuit jpf at time, times_variables.t_pulse_cumulative(3), due to changes in coil currents
                    vpfi[jpf] = vpfi[jpf] + vpfij

                    #  MVA in circuit jpf at time, times_variables.t_pulse_cumulative(3) due to changes in current
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
                # 'time' is the time INDEX.  't_pulse_cumulative' is the time.
                for time in range(6):
                    poloidalenergy[time] = (
                        poloidalenergy[time]
                        + 0.5e0
                        * inductxcurrent[time]
                        * pfcoil_variables.c_pf_coil_turn[jpf, time]
                    )

                #   do time = 1,5
                #     # Mean rate of change of stored energy between time and time+1
                #     if(abs(times_variables.t_pulse_cumulative(time+1)-times_variables.t_pulse_cumulative(time)).gt.1.0e0) :
                #         pf_power_variables.poloidalpower(time) = (poloidalenergy(time+1)-poloidalenergy(time)) / (times_variables.t_pulse_cumulative(time+1)-times_variables.t_pulse_cumulative(time))
                #     else:
                #         # Flag when an interval is small or zero MDK 30/11/16
                #         pf_power_variables.poloidalpower(time) = 9.9e9
                #

                #   end do
                #   #engxpc = 0.5e0 * engx * pfcoil_variables.c_pf_coil_turn(jpf,5)
                #   #ensxpf = ensxpf + engxpc

                #  Resistive power in circuits at times times_variables.t_pulse_cumulative(3) and times_variables.t_pulse_cumulative(5) respectively (MW)
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
            # 'time' is the time INDEX.  't_pulse_cumulative' is the time.
            # Mean rate of change of stored energy between time and time+1
            if (
                abs(
                    times_variables.t_pulse_cumulative[time + 1]
                    - times_variables.t_pulse_cumulative[time]
                )
                > 1.0e0
            ):
                pf_power_variables.poloidalpower[time] = (
                    poloidalenergy[time + 1] - poloidalenergy[time]
                ) / (
                    times_variables.t_pulse_cumulative[time + 1]
                    - times_variables.t_pulse_cumulative[time]
                )
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
        if (
            times_variables.t_pulse_cumulative[4]
            - times_variables.t_pulse_cumulative[3]
            > 1.0e0
        ):
            pfpower = sum(pfdissipation[:]) / (
                times_variables.t_pulse_cumulative[4]
                - times_variables.t_pulse_cumulative[3]
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
        # p_pf_electric_supplies_mw = physics_variables.p_plasma_ohmic_mw / pfcoil_variables.etapsu
        wall_plug_ohmicmw = physics_variables.p_plasma_ohmic_mw * (
            1.0e0 / pfcoil_variables.etapsu - 1.0e0
        )
        # Total mean wall plug power dissipated in PFC and CS power supplies.  Issue #713
        pfcoil_variables.p_pf_electric_supplies_mw = wall_plug_ohmicmw + pfpowermw

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

        # write(self.outfile,50)(times_variables.t_pulse_cumulative(time),time=1,6)

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
        ptfmw = heat_transport_variables.p_tf_electric_supplies_mw
        ppfmw = 1.0e-3 * pf_power_variables.srcktpm
        if pf_power_variables.i_pf_energy_storage_source == 2:
            ppfmw = ppfmw + heat_transport_variables.peakmva

        #  Power to plasma heating supplies, MW
        pheatingmw = (
            heat_transport_variables.p_hcd_electric_total_mw
        )  # Should be zero if i_plasma_ignited==1

        #  Power to cryogenic comp. motors, MW
        crymw = heat_transport_variables.p_cryo_plant_electric_mw

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
        if pf_power_variables.i_pf_energy_storage_source != 2:
            heat_transport_variables.pacpmw = (
                heat_transport_variables.pacpmw + heat_transport_variables.fmgdmw
            )

        # Estimate of the total low voltage power, MW
        # MDK No idea what this is - especially the last term
        # It is used in the old cost routine, so I will leave it in place.
        heat_transport_variables.tlvpmw = (
            heat_transport_variables.p_plant_electric_base_total_mw
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
            "(p_plant_electric_base_total_mw)",
            heat_transport_variables.p_plant_electric_base_total_mw,
            "OP ",
        )
        # MDK Remove this output: no idea what this is
        # po.ovarre(self.outfile,'Total low voltage power (MW)','(tlvpmw)',tlvpmw)

    def component_thermal_powers(self):
        """
        Calculates the first part of the heat transport
        and plant power balance constituents
        author: P J Knight, CCFE, Culham Science Centre
        None
        This routine calculates the first part of the heat transport
        and plant power balance constituents.
        None
        """
        if int(fwbs_variables.i_p_coolant_pumping) not in (2, 3):
            primary_pumping_variables.p_fw_blkt_coolant_pump_mw = (
                heat_transport_variables.p_fw_coolant_pump_mw
                + heat_transport_variables.p_blkt_coolant_pump_mw
            )

        #  Account for pump electrical inefficiencies. The coolant pumps are not assumed to be
        #  100% efficient so the electric power to run them is greater than the power deposited
        #  in the coolant.  The difference should be lost as secondary heat.

        power_variables.p_fw_blkt_coolant_pump_elec_mw = (
            primary_pumping_variables.p_fw_blkt_coolant_pump_mw
            / fwbs_variables.eta_coolant_pump_electric
        )
        power_variables.p_shld_coolant_pump_elec_mw = (
            heat_transport_variables.p_shld_coolant_pump_mw
            / fwbs_variables.eta_coolant_pump_electric
        )
        power_variables.p_div_coolant_pump_elec_mw = (
            heat_transport_variables.p_div_coolant_pump_mw
            / fwbs_variables.eta_coolant_pump_electric
        )

        # Secondary breeder coolant loop. Should return zero if not used.
        power_variables.p_blkt_breeder_pump_elec_mw = (
            heat_transport_variables.p_blkt_breeder_pump_mw
            / fwbs_variables.eta_coolant_pump_electric
        )

        # Total mechanical pump power needed (deposited in coolant)
        power_variables.p_coolant_pump_total_mw = (
            primary_pumping_variables.p_fw_blkt_coolant_pump_mw
            + heat_transport_variables.p_blkt_breeder_pump_mw
            + heat_transport_variables.p_shld_coolant_pump_mw
            + heat_transport_variables.p_div_coolant_pump_mw
        )

        # Minimum total electrical power for primary coolant pumps (MW)
        heat_transport_variables.p_coolant_pump_elec_total_mw = (
            power_variables.p_fw_blkt_coolant_pump_elec_mw
            + power_variables.p_blkt_breeder_pump_elec_mw
            + power_variables.p_shld_coolant_pump_elec_mw
            + power_variables.p_div_coolant_pump_elec_mw
        )

        #  Heat lost through pump power inefficiencies (MW)
        heat_transport_variables.p_coolant_pump_loss_total_mw = (
            heat_transport_variables.p_coolant_pump_elec_total_mw
            - power_variables.p_coolant_pump_total_mw
        )

        # Heat lost in power supplies for heating and current drive
        heat_transport_variables.p_hcd_electric_loss_mw = (
            heat_transport_variables.p_hcd_electric_total_mw
            - current_drive_variables.p_hcd_injected_total_mw
        )

        # Liquid metal breeder/coolant
        # Calculate fraction of blanket nuclear power deposited in liquid breeder / coolant
        if fwbs_variables.i_blkt_dual_coolant == 2:
            power_variables.p_blkt_liquid_breeder_heat_deposited_mw = (
                fwbs_variables.p_blkt_nuclear_heat_total_mw
                * fwbs_variables.f_nuc_pow_bz_liq
            ) + heat_transport_variables.p_blkt_breeder_pump_mw

        # Liquid breeder is circulated but does no cooling
        elif fwbs_variables.i_blkt_dual_coolant == 1:
            power_variables.p_blkt_liquid_breeder_heat_deposited_mw = (
                heat_transport_variables.p_blkt_breeder_pump_mw
            )

        # Liquid breeder also acts a coolant
        if int(fwbs_variables.i_blkt_dual_coolant) in [1, 2]:
            power_variables.p_fw_blkt_heat_deposited_mw = (
                fwbs_variables.p_fw_nuclear_heat_total_mw
                + fwbs_variables.p_fw_rad_total_mw
                + fwbs_variables.p_blkt_nuclear_heat_total_mw
                + heat_transport_variables.p_blkt_breeder_pump_mw
                + primary_pumping_variables.p_fw_blkt_coolant_pump_mw
                + current_drive_variables.p_beam_orbit_loss_mw
                + physics_variables.p_fw_alpha_mw
                + current_drive_variables.p_beam_shine_through_mw
            )
        else:
            # No secondary liquid metal breeder/coolant
            power_variables.p_fw_blkt_heat_deposited_mw = (
                fwbs_variables.p_fw_nuclear_heat_total_mw
                + fwbs_variables.p_fw_rad_total_mw
                + fwbs_variables.p_blkt_nuclear_heat_total_mw
                + primary_pumping_variables.p_fw_blkt_coolant_pump_mw
                + current_drive_variables.p_beam_orbit_loss_mw
                + physics_variables.p_fw_alpha_mw
                + current_drive_variables.p_beam_shine_through_mw
            )

        #  Total power deposited in first wall coolant (MW)
        power_variables.p_fw_heat_deposited_mw = (
            fwbs_variables.p_fw_nuclear_heat_total_mw
            + fwbs_variables.p_fw_rad_total_mw
            + heat_transport_variables.p_fw_coolant_pump_mw
            + current_drive_variables.p_beam_orbit_loss_mw
            + physics_variables.p_fw_alpha_mw
            + current_drive_variables.p_beam_shine_through_mw
        )

        #  Total power deposited in blanket coolant (MW)
        power_variables.p_blkt_heat_deposited_mw = (
            fwbs_variables.p_blkt_nuclear_heat_total_mw
            + heat_transport_variables.p_blkt_coolant_pump_mw
        )

        #  Total power deposited in shield coolant (MW)
        power_variables.p_shld_heat_deposited_mw = (
            fwbs_variables.p_cp_shield_nuclear_heat_mw
            + fwbs_variables.p_shld_nuclear_heat_mw
            + heat_transport_variables.p_shld_coolant_pump_mw
        )

        #  Total thermal power deposited in divertor (MW)
        power_variables.p_div_heat_deposited_mw = (
            physics_variables.p_plasma_separatrix_mw
            + (
                fwbs_variables.p_div_nuclear_heat_total_mw
                + fwbs_variables.p_div_rad_total_mw
            )
            + heat_transport_variables.p_div_coolant_pump_mw
        )

        #  Heat removal from first wall and divertor (MW) (only used in costs.f90)
        if fwbs_variables.i_p_coolant_pumping != 3:
            heat_transport_variables.p_fw_div_heat_deposited_mw = (
                power_variables.p_fw_heat_deposited_mw
                + power_variables.p_div_heat_deposited_mw
            )

        #  Thermal to electric efficiency
        heat_transport_variables.eta_turbine = self.plant_thermal_efficiency(
            heat_transport_variables.eta_turbine
        )
        heat_transport_variables.etath_liq = self.plant_thermal_efficiency_2(
            heat_transport_variables.etath_liq
        )

        #  Primary (high-grade) thermal power, available for electricity generation.  Switch heat_transport_variables.i_shld_primary_heat
        #  is 1 or 0, is user choice on whether the shield thermal power goes to primary or secondary heat
        if fwbs_variables.i_thermal_electric_conversion == 0:
            #  Primary thermal power (MW)
            heat_transport_variables.p_plant_primary_heat_mw = (
                power_variables.p_fw_blkt_heat_deposited_mw
                + heat_transport_variables.i_shld_primary_heat
                * power_variables.p_shld_heat_deposited_mw
            )
            #  Secondary thermal power deposited in divertor (MW)
            heat_transport_variables.p_div_secondary_heat_mw = (
                power_variables.p_div_heat_deposited_mw
            )
            # Divertor primary/secondary power switch: does NOT contribute to energy generation cycle
            power_variables.i_div_primary_heat = 0
        else:
            #  Primary thermal power used to generate electricity (MW)
            heat_transport_variables.p_plant_primary_heat_mw = (
                power_variables.p_fw_blkt_heat_deposited_mw
                + heat_transport_variables.i_shld_primary_heat
                * power_variables.p_shld_heat_deposited_mw
                + power_variables.p_div_heat_deposited_mw
            )
            #  Secondary thermal power deposited in divertor (MW)
            heat_transport_variables.p_div_secondary_heat_mw = 0.0e0
            # Divertor primary/secondary power switch: contributes to energy generation cycle
            power_variables.i_div_primary_heat = 1

        if abs(heat_transport_variables.p_plant_primary_heat_mw) < 1.0e-4:
            logger.error(f"{'ERROR Primary thermal power is zero or negative'}")

        # #284 Fraction of total high-grade thermal power to divertor
        power_variables.f_p_div_primary_heat = (
            power_variables.p_div_heat_deposited_mw
            / heat_transport_variables.p_plant_primary_heat_mw
        )
        # Loss in efficiency as this primary power is collecetd at very low temperature
        power_variables.delta_eta = 0.339 * power_variables.f_p_div_primary_heat

        # ===============================================
        #  Secondary thermal powers
        # ================================================

        #  Secondary thermal power deposited in shield
        heat_transport_variables.p_shld_secondary_heat_mw = (
            power_variables.p_shld_heat_deposited_mw
            * (1 - heat_transport_variables.i_shld_primary_heat)
        )

        #  Secondary thermal power lost to HCD apparatus and diagnostics
        heat_transport_variables.p_hcd_secondary_heat_mw = (
            fwbs_variables.p_fw_hcd_nuclear_heat_mw
            + fwbs_variables.p_fw_hcd_rad_total_mw
        )

        #  Number of primary heat exchangers
        heat_transport_variables.n_primary_heat_exchangers = math.ceil(
            heat_transport_variables.p_plant_primary_heat_mw / 1000.0e0
        )

    def calculate_cryo_loads(self) -> None:
        """
        Calculates and updates the cryogenic heat loads for the system.

        This method computes the various cryogenic heat loads, including conduction/radiation,
        nuclear heating, AC losses, and resistive losses in current leads. It also updates
        the miscellaneous allowance and total heat removal at cryogenic temperatures.
        The results are stored in the corresponding instance variables.
        """

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
                fwbs_variables.p_tf_nuclear_heat_mw,
                pf_power_variables.ensxpfm,
                times_variables.t_plant_pulse_plasma_present,
                tfcoil_variables.c_tf_turn,
                tfcoil_variables.n_tf_coils,
            )

            # Use 13% of ideal Carnot efficiency to fit J. Miller estimate
            # Rem SK : This ITER efficiency is very low compare to the Strowbridge curve
            #          any reasons why?
            # Calculate electric power requirement for cryogenic plant at tfcoil_variables.temp_tf_cryo (MW)
            heat_transport_variables.p_cryo_plant_electric_mw = (
                1.0e-6
                * (constants.TEMP_ROOM - tfcoil_variables.temp_tf_cryo)
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
                + tfcoil_variables.p_tf_joints_resistive
                + fwbs_variables.pnuc_cp_tf * 1.0e6
            )

            # Calculate electric power requirement for cryogenic plant at tfcoil_variables.tcoolin (MW)
            p_tf_cryoal_cryo = (
                1.0e-6
                * (constants.TEMP_ROOM - tfcoil_variables.tcoolin)
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

    def output_plant_thermal_powers(self):
        po.oheadr(self.outfile, "Plant Heat Transport Balance")

        po.ocmmnt(self.outfile, "First Wall : ")
        po.oblnkl(self.outfile)
        po.ovarre(
            self.outfile,
            "Neutronic nuclear heat deposited in FW [MW]",
            "(p_fw_nuclear_heat_total_mw)",
            fwbs_variables.p_fw_nuclear_heat_total_mw,
        )
        po.ovarre(
            self.outfile,
            "Radiation heat deposited in FW [MW]",
            "(p_fw_rad_total_mw)",
            fwbs_variables.p_fw_rad_total_mw,
        )
        po.ovarre(
            self.outfile,
            "Lost alpha-particle heat deposited in FW [MW]",
            "(p_fw_alpha_mw)",
            physics_variables.p_fw_alpha_mw,
        )
        po.ovarre(
            self.outfile,
            "Neutral beam shine-through heat deposited in FW [MW]",
            "(p_beam_shine_through_mw)",
            current_drive_variables.p_beam_shine_through_mw,
        )
        po.ovarre(
            self.outfile,
            "Neutral beam orbit loss heat deposited in FW [MW]",
            "(p_beam_orbit_loss_mw)",
            current_drive_variables.p_beam_orbit_loss_mw,
        )

        po.ovarre(
            self.outfile,
            "Mechancial pumping power deposited in FW coolant [MW]",
            "(p_fw_coolant_pump_mw)",
            heat_transport_variables.p_fw_coolant_pump_mw,
        )
        po.oblnkl(self.outfile)
        po.ovarre(
            self.outfile,
            "Total heat deposited in FW and coolant [MW]",
            "(p_fw_heat_deposited_mw)",
            power_variables.p_fw_heat_deposited_mw,
        )
        po.oblnkl(self.outfile)

        po.ocmmnt(self.outfile, "Blanket : ")
        po.oblnkl(self.outfile)
        po.ovarre(
            self.outfile,
            "Total neutronic nuclear heat deposited and created in Blanket(s) [MW]",
            "(p_blkt_nuclear_heat_total_mw)",
            fwbs_variables.p_blkt_nuclear_heat_total_mw,
        )
        po.ovarre(
            self.outfile,
            "Total multiplication neutronic nuclear heat created in Blanket(s) [MW]",
            "(p_blkt_multiplication_mw)",
            fwbs_variables.p_blkt_multiplication_mw,
        )
        po.ovarre(
            self.outfile,
            "Neutron nuclear heat multiplication factor in Blanket(s)",
            "(f_p_blkt_multiplication)",
            fwbs_variables.f_p_blkt_multiplication,
        )

        po.ovarre(
            self.outfile,
            "Mechancial pumping power deposited in Blanket(s) coolant [MW]",
            "(p_blkt_coolant_pump_mw)",
            heat_transport_variables.p_blkt_coolant_pump_mw,
        )

        po.oblnkl(self.outfile)
        po.ovarre(
            self.outfile,
            "Total heat deposited in Blanket(s) and coolant [MW]",
            "(p_blkt_heat_deposited_mw)",
            power_variables.p_blkt_heat_deposited_mw,
        )

        po.oblnkl(self.outfile)
        po.ocmmnt(self.outfile, "FW and Blanket : ")
        po.oblnkl(self.outfile)

        po.ovarre(
            self.outfile,
            "Mechancial pumping power deposited in Blanket(s) and FW coolant [MW]",
            "(p_fw_blkt_coolant_pump_mw)",
            primary_pumping_variables.p_fw_blkt_coolant_pump_mw,
        )
        po.ovarre(
            self.outfile,
            "Total heat deposited in Blanket(s) and FW coolant [MW]",
            "(p_fw_blkt_heat_deposited_mw)",
            power_variables.p_fw_blkt_heat_deposited_mw,
        )

        po.oblnkl(self.outfile)
        po.ocmmnt(self.outfile, "----------------------------")
        po.oblnkl(self.outfile)
        po.ocmmnt(self.outfile, "VV and Shield : ")
        po.oblnkl(self.outfile)
        po.ovarre(
            self.outfile,
            "Neutronic nuclear heat deposited in VV shield [MW]",
            "(p_shld_nuclear_heat_mw)",
            fwbs_variables.p_shld_nuclear_heat_mw,
        )
        po.ovarre(
            self.outfile,
            "Neutronic nuclear heat deposited in ST centrepost shield [MW]",
            "(p_cp_shield_nuclear_heat_mw)",
            fwbs_variables.p_cp_shield_nuclear_heat_mw,
        )

        po.ovarre(
            self.outfile,
            "Mechancial pumping power deposited in shield coolant(s) [MW]",
            "(p_shld_coolant_pump_mw)",
            heat_transport_variables.p_shld_coolant_pump_mw,
        )

        po.oblnkl(self.outfile)
        po.ovarre(
            self.outfile,
            "Total heat deposited in VV and shield coolant(s) [MW]",
            "(p_shld_heat_deposited_mw)",
            power_variables.p_shld_heat_deposited_mw,
        )

        po.oblnkl(self.outfile)
        po.ocmmnt(self.outfile, "----------------------------")
        po.oblnkl(self.outfile)
        po.ocmmnt(self.outfile, "Divertor : ")
        po.oblnkl(self.outfile)
        po.ovarre(
            self.outfile,
            "Plasma separatrix power deposited in divertor [MW]",
            "(p_plasma_separatrix_mw)",
            physics_variables.p_plasma_separatrix_mw,
        )
        po.ovarre(
            self.outfile,
            "Neutronic nuclear heat deposited in divertor [MW]",
            "(p_div_nuclear_heat_total_mw)",
            fwbs_variables.p_div_nuclear_heat_total_mw,
        )
        po.ovarre(
            self.outfile,
            "Radiation heat deposited in divertor [MW]",
            "(p_div_rad_total_mw)",
            fwbs_variables.p_div_rad_total_mw,
        )

        po.ovarre(
            self.outfile,
            "Mechancial pumping power deposited in divertor coolant [MW]",
            "(p_div_coolant_pump_mw)",
            heat_transport_variables.p_div_coolant_pump_mw,
        )

        po.oblnkl(self.outfile)
        po.ovarre(
            self.outfile,
            "Total heat deposited in divertor and coolants [MW]",
            "(p_div_heat_deposited_mw)",
            power_variables.p_div_heat_deposited_mw,
        )

        po.oblnkl(self.outfile)
        po.ocmmnt(self.outfile, "----------------------------")
        po.oblnkl(self.outfile)

        po.ovarre(
            self.outfile,
            "Mechanical pumping power of all coolant pumps [MW]",
            "(p_coolant_pump_total_mw)",
            power_variables.p_coolant_pump_total_mw,
        )

        po.oblnkl(self.outfile)
        po.ocmmnt(self.outfile, "----------------------------")
        po.oblnkl(self.outfile)
        po.ocmmnt(self.outfile, "Secondary heat : ")
        po.oblnkl(self.outfile)

        po.ovarre(
            self.outfile,
            "Electric power for core plant systems [MW]",
            "(p_plant_core_systems_elec_mw)",
            power_variables.p_plant_core_systems_elec_mw,
        )
        po.ovarre(
            self.outfile,
            "Wall plug losses in H&CD systems [MW]",
            "(p_hcd_electric_loss_mw)",
            heat_transport_variables.p_hcd_electric_loss_mw,
        )
        po.ovarre(
            self.outfile,
            "Total wall plug losses in coolant pump systems [MW]",
            "(p_coolant_pump_loss_total_mw)",
            heat_transport_variables.p_coolant_pump_loss_total_mw,
        )
        po.ovarre(
            self.outfile,
            "Divertor thermal power not used for electricity production [MW]",
            "(p_div_secondary_heat_mw)",
            heat_transport_variables.p_div_secondary_heat_mw,
        )
        po.ovarre(
            self.outfile,
            "Shield thermal power not used for electricity production [MW]",
            "(p_shld_secondary_heat_mw)",
            heat_transport_variables.p_shld_secondary_heat_mw,
        )
        po.ovarre(
            self.outfile,
            "Neutronic nuclear heating in TF coils [MW]",
            "(p_tf_nuclear_heat_mw)",
            fwbs_variables.p_tf_nuclear_heat_mw,
        )
        po.oblnkl(self.outfile)
        po.ovarre(
            self.outfile,
            "Neutronic nuclear heating in H&CD systems and diagnostics [MW]",
            "(p_fw_hcd_nuclear_heat_mw)",
            fwbs_variables.p_fw_hcd_nuclear_heat_mw,
        )
        po.ovarre(
            self.outfile,
            "Radiation heat deposited in H&CD systems and diagnostics [MW]",
            "(p_fw_hcd_rad_total_mw)",
            fwbs_variables.p_fw_hcd_rad_total_mw,
        )
        po.ovarre(
            self.outfile,
            "Total heat deposited in in H&CD systems and diagnostics [MW]",
            "(p_hcd_secondary_heat_mw)",
            heat_transport_variables.p_hcd_secondary_heat_mw,
        )
        po.oblnkl(self.outfile)
        po.ovarre(
            self.outfile,
            "Total secondary heat not used for electricity production [MW]",
            "(p_plant_secondary_heat_mw)",
            heat_transport_variables.p_plant_secondary_heat_mw,
        )

        po.oblnkl(self.outfile)
        po.ocmmnt(self.outfile, "----------------------------")
        po.oblnkl(self.outfile)
        po.ocmmnt(self.outfile, "Primary heat : ")
        po.oblnkl(self.outfile)

        po.ovarre(
            self.outfile,
            "Total heat deposited in FW and coolant [MW]",
            "(p_fw_heat_deposited_mw)",
            power_variables.p_fw_heat_deposited_mw,
        )
        po.ovarre(
            self.outfile,
            "Total heat deposited in Blanket(s) and coolant [MW]",
            "(p_blkt_heat_deposited_mw)",
            power_variables.p_blkt_heat_deposited_mw,
        )
        po.ovarre(
            self.outfile,
            "Total heat deposited in Blanket(s) and FW coolant [MW]",
            "(p_fw_blkt_heat_deposited_mw)",
            power_variables.p_fw_blkt_heat_deposited_mw,
        )
        po.oblnkl(self.outfile)

        po.ovarre(
            self.outfile,
            "Total heat deposited in VV and shield coolant(s) [MW]",
            "(p_shld_heat_deposited_mw)",
            power_variables.p_shld_heat_deposited_mw,
        )
        po.oblnkl(self.outfile)
        po.ovarre(
            self.outfile,
            "Total heat deposited in divertor and coolants [MW]",
            "(p_div_heat_deposited_mw)",
            power_variables.p_div_heat_deposited_mw,
        )
        po.ovarre(
            self.outfile,
            "Fraction of total primary heat originating from divertor",
            "(f_p_div_primary_heat)",
            power_variables.f_p_div_primary_heat,
        )
        po.oblnkl(self.outfile)
        po.ovarre(
            self.outfile,
            "Total primary thermal power used for electricity production [MW]",
            "(p_plant_primary_heat_mw)",
            heat_transport_variables.p_plant_primary_heat_mw,
        )

    def output_plant_electric_powers(self):
        po.oheadr(self.outfile, "Plant Electricity Production")

        po.ocmmnt(self.outfile, "Turbine conversion : ")
        po.oblnkl(self.outfile)

        po.ovarre(
            self.outfile,
            "Total high grade thermal power used for electricity production [MWth]",
            "(p_plant_primary_heat_mw)",
            heat_transport_variables.p_plant_primary_heat_mw,
        )

        po.ovarrf(
            self.outfile,
            "Thermal to electric conversion efficiency of the turbine",
            "(eta_turbine)",
            heat_transport_variables.eta_turbine,
        )
        po.ovarre(
            self.outfile,
            "Total thermal power lost in power conversion [MWth]",
            "(p_turbine_loss_mw)",
            power_variables.p_turbine_loss_mw,
        )
        po.oblnkl(self.outfile)
        po.ovarre(
            self.outfile,
            "Total electric power produced [MWe]",
            "(p_plant_electric_gross_mw)",
            heat_transport_variables.p_plant_electric_gross_mw,
        )
        po.oblnkl(self.outfile)
        po.ocmmnt(self.outfile, "----------------------------")
        po.oblnkl(self.outfile)

        po.ocmmnt(self.outfile, "Electric requirements of core plant systems : ")
        po.oblnkl(self.outfile)

        po.ovarre(
            self.outfile,
            "Base plant electric load [We]",
            "(p_plant_electric_base)",
            heat_transport_variables.p_plant_electric_base,
        )
        po.ovarre(
            self.outfile,
            "Electric power per unit area of plant floor space [We/m^2]",
            "(pflux_plant_floor_electric)",
            heat_transport_variables.pflux_plant_floor_electric,
        )
        po.ovarre(
            self.outfile,
            "Effective area of plant buildings floor [m^2]",
            "(a_plant_floor_effective)",
            buildings_variables.a_plant_floor_effective,
        )
        po.oblnkl(self.outfile)
        po.ovarre(
            self.outfile,
            "Total base plant electric load [MWe]",
            "(p_plant_electric_base_total_mw)",
            heat_transport_variables.p_plant_electric_base_total_mw,
        )
        po.oblnkl(self.outfile)
        po.ovarre(
            self.outfile,
            "Electric power demand for cryo plant [MWe]",
            "(p_cryo_plant_electric_mw)",
            heat_transport_variables.p_cryo_plant_electric_mw,
        )
        po.ovarre(
            self.outfile,
            "Electric power demand for tritium plant [MWe]",
            "(p_tritium_plant_electric_mw)",
            heat_transport_variables.p_tritium_plant_electric_mw,
        )
        po.ovarre(
            self.outfile,
            "Electric power demand for vacuum pumps [MWe]",
            "(vachtmw)",
            heat_transport_variables.vachtmw,
        )
        po.ovarre(
            self.outfile,
            "Electric power demand for TF coil system [MWe]",
            "(p_tf_electric_supplies_mw)",
            heat_transport_variables.p_tf_electric_supplies_mw,
        )
        po.ovarre(
            self.outfile,
            "Electric power demand for PF coil system [MWe]",
            "(p_pf_electric_supplies_mw)",
            pfcoil_variables.p_pf_electric_supplies_mw,
        )
        po.ovarre(
            self.outfile,
            "Electric power demand for CP coolant pumps [MWe]",
            "(p_cp_coolant_pump_elec_mw)",
            power_variables.p_cp_coolant_pump_elec_mw,
        )
        po.oblnkl(self.outfile)
        po.ovarre(
            self.outfile,
            "Electric power demand of core plant systems needed at all times [MWe]",
            "(p_plant_core_systems_elec_mw)",
            power_variables.p_plant_core_systems_elec_mw,
        )

        po.oblnkl(self.outfile)
        po.ocmmnt(self.outfile, "----------------------------")
        po.oblnkl(self.outfile)

        po.ocmmnt(self.outfile, "Electric requirements during plasma flat-top : ")
        po.oblnkl(self.outfile)

        po.ovarre(
            self.outfile,
            "Electric power demand of FW and Blanket coolant pumps [MWe]",
            "(p_fw_blkt_coolant_pump_elec_mw)",
            power_variables.p_fw_blkt_coolant_pump_elec_mw,
        )
        po.ovarre(
            self.outfile,
            "Electric power demand of Blanket secondary breeder coolant pumps [MWe]",
            "(p_blkt_breeder_pump_elec_mw)",
            power_variables.p_blkt_breeder_pump_elec_mw,
        )
        po.ovarre(
            self.outfile,
            "Electric power demand of VV and Shield coolant pumps [MWe]",
            "(p_shld_coolant_pump_elec_mw)",
            power_variables.p_shld_coolant_pump_elec_mw,
        )
        po.ovarre(
            self.outfile,
            "Electric power demand of Divertor colant pumps [MWe]",
            "(p_div_coolant_pump_elec_mw)",
            power_variables.p_div_coolant_pump_elec_mw,
        )
        po.oblnkl(self.outfile)
        po.ovarre(
            self.outfile,
            "Electric wall plug efficiency of coolant pumps",
            "(eta_coolant_pump_electric)",
            fwbs_variables.eta_coolant_pump_electric,
        )
        po.ovarre(
            self.outfile,
            "Total electric demand of all coolant pumps [MWe]",
            "(p_coolant_pump_elec_total_mw)",
            heat_transport_variables.p_coolant_pump_elec_total_mw,
        )
        po.oblnkl(self.outfile)
        po.ovarre(
            self.outfile,
            "Total electric demand of all H&CD systems [MWe]",
            "(p_hcd_electric_total_mw)",
            heat_transport_variables.p_hcd_electric_total_mw,
        )

        po.oblnkl(self.outfile)

        po.ovarre(
            self.outfile,
            "Total re-circulated electric power of the plant [MWe]",
            "(p_plant_electric_recirc_mw)",
            heat_transport_variables.p_plant_electric_recirc_mw,
        )
        po.ovarre(
            self.outfile,
            "Fraction of gross electricity re-circulated",
            "(f_p_plant_electric_recirc)",
            heat_transport_variables.f_p_plant_electric_recirc,
        )

        po.oblnkl(self.outfile)

        po.ovarre(
            self.outfile,
            "Total net-electric power of the plant [MWe]",
            "(p_plant_electric_net_mw)",
            heat_transport_variables.p_plant_electric_net_mw,
        )

        po.oblnkl(self.outfile)

        po.ovarre(
            self.outfile,
            "Total electric energy output per pulse (MJ)",
            "(e_plant_net_electric_pulse_mj)",
            power_variables.e_plant_net_electric_pulse_mj,
        )
        po.ovarre(
            self.outfile,
            "Total electric energy output per pulse (kWh)",
            "(e_plant_net_electric_pulse_kwh)",
            power_variables.e_plant_net_electric_pulse_kwh,
        )

    def plant_electric_production(self) -> None:
        """
        This method completes the calculation of the plant's electrical and thermal power flows,
        including secondary heat, recirculating power, net and gross electric power, and various
        efficiency measures.

        If `output` is True, the method writes a comprehensive summary of the plant's power and
        heat transport balance, assumptions, and efficiency metrics to the specified output file.

        """
        if physics_variables.itart == 1 and tfcoil_variables.i_tf_sup == 0:
            power_variables.p_cp_coolant_pump_elec_mw = (
                1.0e-6 * tfcoil_variables.p_cp_coolant_pump_elec
            )
        else:
            power_variables.p_cp_coolant_pump_elec_mw = 0.0e0

        #  Total baseline power to facility loads, MW
        heat_transport_variables.p_plant_electric_base_total_mw = (
            heat_transport_variables.p_plant_electric_base * 1.0e-6
            + buildings_variables.a_plant_floor_effective
            * (heat_transport_variables.pflux_plant_floor_electric * 1.0e-3)
            / 1000.0e0
        )

        #  Facility heat removal (heat_transport_variables.p_plant_electric_base_total_mw calculated in ACPOW)
        heat_transport_variables.fachtmw = (
            heat_transport_variables.p_plant_electric_base_total_mw
        )

        #  Electrical power consumed by fusion power core systems
        #  (excluding heat transport pumps and auxiliary injection power system)

        power_variables.p_plant_core_systems_elec_mw = (
            heat_transport_variables.p_cryo_plant_electric_mw
            + heat_transport_variables.fachtmw
            + power_variables.p_cp_coolant_pump_elec_mw
            + heat_transport_variables.p_tf_electric_supplies_mw
            + heat_transport_variables.p_tritium_plant_electric_mw
            + heat_transport_variables.vachtmw
            + pfcoil_variables.p_pf_electric_supplies_mw
        )

        #  Total secondary heat
        #  (total low-grade heat rejected - does not contribute to power conversion cycle)
        #  Included fwbs_variables.p_tf_nuclear_heat_mw
        # p_plant_secondary_heat_mw = power_variables.p_plant_core_systems_elec_mw + heat_transport_variables.p_hcd_electric_loss_mw + heat_transport_variables.p_coolant_pump_loss_total_mw + hthermmw + heat_transport_variables.p_div_secondary_heat_mw + heat_transport_variables.p_shld_secondary_heat_mw + heat_transport_variables.p_hcd_secondary_heat_mw + fwbs_variables.p_tf_nuclear_heat_mw
        heat_transport_variables.p_plant_secondary_heat_mw = (
            power_variables.p_plant_core_systems_elec_mw
            + heat_transport_variables.p_hcd_electric_loss_mw
            + heat_transport_variables.p_coolant_pump_loss_total_mw
            + heat_transport_variables.p_div_secondary_heat_mw
            + heat_transport_variables.p_shld_secondary_heat_mw
            + heat_transport_variables.p_hcd_secondary_heat_mw
            + fwbs_variables.p_tf_nuclear_heat_mw
        )

        #  Calculate powers relevant to a power-producing plant
        if cost_variables.ireactor == 1:
            #  Gross electric power
            # p_plant_electric_gross_mw = (heat_transport_variables.p_plant_primary_heat_mw-hthermmw) * heat_transport_variables.eta_turbine
            if (
                fwbs_variables.i_blkt_dual_coolant > 0
                and fwbs_variables.i_p_coolant_pumping == 2
            ):
                heat_transport_variables.p_plant_electric_gross_mw = (
                    (
                        heat_transport_variables.p_plant_primary_heat_mw
                        - power_variables.p_blkt_liquid_breeder_heat_deposited_mw
                    )
                    * heat_transport_variables.eta_turbine
                    + power_variables.p_blkt_liquid_breeder_heat_deposited_mw
                    * heat_transport_variables.etath_liq
                )
            else:
                heat_transport_variables.p_plant_electric_gross_mw = (
                    heat_transport_variables.p_plant_primary_heat_mw
                    * heat_transport_variables.eta_turbine
                )

            # Total lost thermal power in the turbine
            power_variables.p_turbine_loss_mw = (
                heat_transport_variables.p_plant_primary_heat_mw
                * (1 - heat_transport_variables.eta_turbine)
            )

            #  Total recirculating power
            heat_transport_variables.p_plant_electric_recirc_mw = (
                power_variables.p_plant_core_systems_elec_mw
                + heat_transport_variables.p_hcd_electric_total_mw
                + heat_transport_variables.p_coolant_pump_elec_total_mw
            )

            #  Net electric power
            heat_transport_variables.p_plant_electric_net_mw = (
                heat_transport_variables.p_plant_electric_gross_mw
                - heat_transport_variables.p_plant_electric_recirc_mw
            )

            #  Recirculating power fraction
            heat_transport_variables.f_p_plant_electric_recirc = (
                heat_transport_variables.p_plant_electric_gross_mw
                - heat_transport_variables.p_plant_electric_net_mw
            ) / heat_transport_variables.p_plant_electric_gross_mw

        (
            power_variables.e_plant_net_electric_pulse_kwh,
            power_variables.e_plant_net_electric_pulse_mj,
            power_variables.p_plant_electric_base_total_profile_mw,
            power_variables.p_plant_electric_gross_profile_mw,
            power_variables.p_plant_electric_net_profile_mw,
            power_variables.p_hcd_electric_total_profile_mw,
            power_variables.p_coolant_pump_elec_total_profile_mw,
            power_variables.p_tf_electric_supplies_profile_mw,
            power_variables.p_pf_electric_supplies_profile_mw,
            power_variables.vachtmw_profile_mw,
            power_variables.p_tritium_plant_electric_profile_mw,
            power_variables.p_cryo_plant_electric_profile_mw,
            power_variables.p_fusion_total_profile_mw,
        ) = self.power_profiles_over_time(
            t_precharge=times_variables.t_plant_pulse_coil_precharge,
            t_current_ramp_up=times_variables.t_plant_pulse_plasma_current_ramp_up,
            t_fusion_ramp=times_variables.t_plant_pulse_fusion_ramp,
            t_burn=times_variables.t_plant_pulse_burn,
            t_ramp_down=times_variables.t_plant_pulse_plasma_current_ramp_down,
            t_between_pulse=times_variables.t_plant_pulse_dwell,
            p_plant_electric_base_total_mw=heat_transport_variables.p_plant_electric_base_total_mw,
            p_cryo_plant_electric_mw=heat_transport_variables.p_cryo_plant_electric_mw,
            p_tritium_plant_electric_mw=heat_transport_variables.p_tritium_plant_electric_mw,
            vachtmw=heat_transport_variables.vachtmw,
            p_tf_electric_supplies_mw=heat_transport_variables.p_tf_electric_supplies_mw,
            p_pf_electric_supplies_mw=pfcoil_variables.p_pf_electric_supplies_mw,
            p_coolant_pump_elec_total_mw=heat_transport_variables.p_coolant_pump_elec_total_mw,
            p_hcd_electric_total_mw=heat_transport_variables.p_hcd_electric_total_mw,
            p_fusion_total_mw=physics_variables.p_fusion_total_mw,
            p_plant_electric_gross_mw=heat_transport_variables.p_plant_electric_gross_mw,
            p_plant_electric_net_mw=heat_transport_variables.p_plant_electric_net_mw,
        )

    def cryo(
        self,
        i_tf_sup,
        tfcryoarea,
        coldmass,
        p_tf_nuclear_heat_mw,
        ensxpfm,
        t_plant_pulse_plasma_present,
        c_tf_turn,
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
        p_tf_nuclear_heat_mw : input real : Nuclear heating in TF coils (MW)
        ensxpfm : input real : Maximum PF coil stored energy (MJ)
        t_plant_pulse_plasma_present : input real : Pulse length of cycle (s)
        c_tf_turn : input real : Current per turn in TF coils (A)
        tfno : input real : Number of TF coils
        helpow : output real : Helium heat removal at cryo temperatures (W)
        This routine calculates the cryogenic heat load.
        D. Slack memo SCMDG 88-5-1-059, LLNL ITER-88-054, Aug. 1988
        """
        power_variables.qss = 4.3e-4 * coldmass
        if i_tf_sup == 1:
            power_variables.qss = power_variables.qss + 2.0e0 * tfcryoarea

        #  Nuclear heating of TF coils (W) (zero if resistive)
        if fwbs_variables.inuclear == 0 and i_tf_sup == 1:
            fwbs_variables.qnuc = 1.0e6 * p_tf_nuclear_heat_mw
        # Issue #511: if fwbs_variables.inuclear = 1 : fwbs_variables.qnuc is input.

        #  AC losses
        power_variables.qac = 1.0e3 * ensxpfm / t_plant_pulse_plasma_present

        #  Current leads
        if i_tf_sup == 1:
            power_variables.qcl = 13.6e-3 * n_tf_coils * c_tf_turn
        else:
            power_variables.qcl = 0.0e0

        #  45% extra miscellaneous, piping and reserves
        power_variables.qmisc = 0.45e0 * (
            power_variables.qss
            + fwbs_variables.qnuc
            + power_variables.qac
            + power_variables.qcl
        )
        return max(
            0.0e0,
            power_variables.qmisc
            + power_variables.qss
            + fwbs_variables.qnuc
            + power_variables.qac
            + power_variables.qcl,
        )

    def output_cryogenics(self) -> None:
        """
        Outputs cryogenic system heat loads and related parameters to the output file.

        This method prints the breakdown of cryogenic heat loads, including conduction/radiation,
        nuclear heating, AC losses, resistive losses in current leads, miscellaneous allowances,
        and total heat removal at cryogenic temperatures. It also outputs the temperatures and
        efficiencies of the cryogenic systems, as well as the electric power required for the
        cryogenic plant.
        """

        po.oheadr(self.outfile, "Cryogenics")
        po.ovarre(
            self.outfile,
            "Conduction and radiation heat loads on cryogenic components (MW)",
            "(qss/1.0d6)",
            power_variables.qss / 1.0e6,
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
            power_variables.qac / 1.0e6,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Resistive losses in current leads (MW)",
            "(qcl/1.0d6)",
            power_variables.qcl / 1.0e6,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "45% allowance for heat loads in transfer lines, storage tanks etc (MW)",
            "(qmisc/1.0d6)",
            power_variables.qmisc / 1.0e6,
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
        po.ovarre(
            self.outfile,
            "Electric power for cryogenic plant (MW)",
            "(p_cryo_plant_electric_mw)",
            heat_transport_variables.p_cryo_plant_electric_mw,
            "OP ",
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
                eta_turbine = 0.411e0 - power_variables.delta_eta
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
                    logger.warning(
                        "Turbine temperature temp_turbine_coolant_in out of range of validity"
                        f"{heat_transport_variables.temp_turbine_coolant_in=}"
                    )

                eta_turbine = (
                    0.1802e0 * np.log(heat_transport_variables.temp_turbine_coolant_in)
                    - 0.7823
                    - power_variables.delta_eta
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
                logger.warning(
                    "Turbine temperature temp_turbine_coolant_in out of range of validity"
                    f"{heat_transport_variables.temp_turbine_coolant_in=}"
                )

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
                logger.warning(
                    "Turbine temperature temp_turbine_coolant_in out of range of validity"
                    f"{heat_transport_variables.temp_turbine_coolant_in=}"
                )

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
            # tfcoil_variables.c_tf_turn  - current per TFC turn (A)
            # tfcoil_variables.j_tf_bus   - bus current density (A/m2)
            a_tf_bus = tfcoil_variables.c_tf_turn / tfcoil_variables.j_tf_bus

            # Bus resistance [ohm]
            # Bus resistivity (tfcoil_variables.rho_tf_bus)
            # Issue #1253: there was a fudge here to set the bus bar resistivity equal
            # to the TF conductor resistivity. I have removed this.
            tfbusres = (
                tfcoil_variables.rho_tf_bus * tfcoil_variables.len_tf_bus / a_tf_bus
            )

            #  Bus mass (kg)
            tfcoil_variables.m_tf_bus = (
                tfcoil_variables.len_tf_bus * a_tf_bus * constants.den_copper
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
                * tfcoil_variables.c_tf_turn
                / tfcoil_variables.n_tf_coils
            )

            # Resistive powers (MW):
            tfcoil_variables.p_cp_resistive_mw = (
                1.0e-6 * tfcoil_variables.p_cp_resistive
            )  # inboard legs (called centrepost, CP for tart design)
            tfcoil_variables.p_tf_leg_resistive_mw = (
                1.0e-6 * tfcoil_variables.p_tf_leg_resistive
            )  # outboard legs
            tfcoil_variables.p_tf_joints_resistive_mw = (
                1.0e-6 * tfcoil_variables.p_tf_joints_resistive
            )  # Joints
            tfbusmw = (
                1.0e-6 * tfcoil_variables.c_tf_turn**2 * tfbusres
            )  # TF coil bus => Dodgy #

            #  TF coil reactive power
            #  Set reactive power to 0, since ramp up can be long
            #  The TF coil can be ramped up as slowly as you like
            #  (although this will affect the time to recover from a magnet quench).
            #     tfreacmw = 1.0e-6 * 1.0e9 * estotf/(t_plant_pulse_plasma_current_ramp_up + t_plant_pulse_coil_precharge)
            #                                 estotf(=e_tf_magnetic_stored_total_gj/tfcoil_variables.n_tf_coils) has been removed (#199 #847)
            tfreacmw = 0.0e0

            # Total power consumption (MW)
            tfcoil_variables.tfcmw = (
                tfcoil_variables.p_cp_resistive_mw
                + tfcoil_variables.p_tf_leg_resistive_mw
                + tfbusmw
                + tfreacmw
                + tfcoil_variables.p_tf_joints_resistive_mw
            )

            # Total steady state AC power demand (MW)
            heat_transport_variables.p_tf_electric_supplies_mw = (
                tfcoil_variables.tfcmw / heat_transport_variables.etatf
            )

        else:  # Superconducting TF coil option
            (
                tfcoil_variables.tfckw,
                tfcoil_variables.len_tf_bus,
                tfcoil_variables.drarea,
                buildings_variables.tfcbv,
                heat_transport_variables.p_tf_electric_supplies_mw,
            ) = self.tfcpwr(
                output=output,
                itfka=tfcoil_variables.c_tf_turn / 1e3,
                rmajor=physics_variables.rmajor,
                n_tf_coils=tfcoil_variables.n_tf_coils,
                v_tf_coil_dump_quench_kv=tfcoil_variables.v_tf_coil_dump_quench_kv,
                e_tf_coil_magnetic_stored_mj=tfcoil_variables.e_tf_coil_magnetic_stored
                / 1e6,
                rptfc=tfcoil_variables.res_tf_leg,
            )
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
            "(p_cp_resistive_mw)",
            tfcoil_variables.p_cp_resistive_mw,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Power dissipation in TF coil set: outboard legs (MW)",
            "(p_tf_leg_resistive_mw)",
            tfcoil_variables.p_tf_leg_resistive_mw,
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
                "(p_tf_joints_resistive_mw)",
                tfcoil_variables.p_tf_joints_resistive_mw,
                "OP ",
            )

        # Reactive poower has been set to zero.
        # po.ovarre(outfile,'TF coil reactive power (MW)','(tfreacmw)', tfreacmw)

    def tfcpwr(
        self,
        output: bool,
        itfka: float,
        rmajor: float,
        n_tf_coils: int,
        v_tf_coil_dump_quench_kv: float,
        e_tf_coil_magnetic_stored_mj: float,
        rptfc: float,
    ) -> tuple[float, float, float, float, float]:
        """
        Calculates the TF coil power conversion system parameters for superconducting coils.

        :param output: If True, outputs results to the output file.
        :type output: bool
        :param itfka: TF coil current (kA).
        :type itfka: float
        :param rmajor: Major radius of the device (m).
        :type rmajor: float
        :param n_tf_coils: Number of TF coils.
        :type n_tf_coils: int
        :param v_tf_coil_dump_quench_kv: Voltage across a TF coil during quench (kV).
        :type v_tf_coil_dump_quench_kv: float
        :param e_tf_coil_magnetic_stored_mj: Stored energy per TF coil (MJ).
        :type e_tf_coil_magnetic_stored_mj: float
        :param rptfc: Resistance per TF coil (ohm).
        :type rptfc: float

        :returns: Tuple containing:
            - tfckw (float): DC power supply rating (kW)
            - len_tf_bus (float): Total length of TF coil bussing (m)
            - drarea (float): Dump resistor floor area (m2)
            - tfcbv (float): TF coil power conversion building volume (m3)
            - p_tf_electric_supplies_mw (float): Total steady state AC power demand (MW)

        :notes:
            - This routine calculates the TF power conversion system parameters: floor space, power supplies,
            bussing, coil protection equipment, and the associated controls and instrumentation.
            Originally written by G. Gorker, FEDC/ORNL, April 1987, modified by J. Galambos in 1991 to run in TETRA,
            and included in PROCESS in 1992 by P. C. Shipe.

            - This routine was originally called tfcpwr() in the ETR/ITER Systems Code

        :references:
            - R. L. Reid, “ETR/ITER Systems Code,” Oak Ridge National Laboratory, ORNL-FEDC-87-7, April 1988.
            Available: https://www.google.com/url?sa=t&rct=j&q=&esrc=s&source=web&cd=&ved=2ahUKEwjX9Ozasb6OAxX1U0EAHV-3C
            0QQFnoECBYQAQ&url=https%3A%2F%2Fengineering.purdue.edu%2FCMUXE%2FPublications%2FAHR%2FR88ORNL-FEDC-87-7.pdf&
            usg=AOvVaw1-LdCefwqI0hJumpHvfTlX&opi=89978449
        ‌
        """

        N_TF_COIL_BREAKERS = 1.0e0  # number of TF coils per circuit breaker
        j_tf_bus_design_ka_cm = 0.125e0  # design current density of TF bus, kA/cm2
        rtfps = 1.05e0  # rating factor for TF coil power supplies
        fspc1 = 0.15e0  # floor space coefficient for power supplies
        fspc2 = 0.8e0  # floor space coefficient for circuit breakers
        fspc3 = 0.4e0  # floor space coefficient for load centres

        if rptfc == 0.0e0:
            T_TF_CHARGE_HOURS = 4.0e0  # charge time of the coils, hours
            nsptfc = 1.0e0  # superconducting (1.0 = superconducting, 0.0 = resistive)
        else:
            T_TF_CHARGE_HOURS = 0.16667e0  # charge time of the coils, hours
            nsptfc = 0.0e0  # resistive (1.0 = superconducting, 0.0 = resistive)

        #  Total steady state TF coil AC power demand (summed later)
        p_tf_electric_supplies_mw = 0.0e0

        #  Stored energy of all TF coils, MJ
        e_tf_magnetic_stored_mj = n_tf_coils * e_tf_coil_magnetic_stored_mj

        #  Inductance of all TF coils, Henries
        ind_tf_total = 2.0e0 * e_tf_magnetic_stored_mj / itfka**2

        #  Number of circuit breakers
        n_tf_breakers = n_tf_coils / N_TF_COIL_BREAKERS

        #  Inductance per TF coil, Henries
        ind_tf_coil = ind_tf_total / n_tf_coils

        #  Aluminium bus section area, sq cm
        albusa = itfka / j_tf_bus_design_ka_cm

        #  Total TF system bus length, m
        len_tf_bus = (
            8.0e0 * np.pi * rmajor
            + (1.0e0 + n_tf_breakers) * (12.0e0 * rmajor + 80.0e0)
            + 0.2e0 * itfka * np.sqrt(n_tf_coils * rptfc * 1000.0e0)
        )

        #  Aluminium bus weight, tonnes
        m_tf_bus_aluminium_tonnes = 2.7e0 * albusa * len_tf_bus / 1.0e4

        #  Total resistance of TF bus, ohms
        # res_tf_bus = 2.62e-4 * len_tf_bus / albusa
        res_tf_bus = tfcoil_variables.rho_tf_bus * len_tf_bus / (albusa / 10000)

        #  Total voltage drop across TF bus, volts
        v_tf_bus = 1000.0e0 * itfka * res_tf_bus

        #  Total resistance of the TF coils, ohms
        rcoils = n_tf_coils * rptfc

        #  Total impedance, ohms
        ztotal = res_tf_bus + rcoils + ind_tf_total / (3600.0e0 * T_TF_CHARGE_HOURS)

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
        r1dump = nsptfc * v_tf_coil_dump_quench_kv * N_TF_COIL_BREAKERS / itfka

        #  Time constant, s
        ttfsec = (
            ind_tf_coil
            * N_TF_COIL_BREAKERS
            / (r1dump * nsptfc + rptfc * (1.0e0 - nsptfc))
        )

        #  Number of dump resistors
        n_tf_dump_resistors = n_tf_breakers * 4.0e0

        #  Peak power to a dump resistor during quench, MW
        r1ppmw = nsptfc * r1dump * (itfka / 2.0e0) ** 2

        #  Energy to dump resistor during quench, MJ
        r1emj = nsptfc * e_tf_magnetic_stored_mj / (n_tf_dump_resistors + 0.0001e0)

        #  Total TF coil peak resistive power demand, MVA
        rpower = (n_tf_coils * rptfc + res_tf_bus) * itfka**2

        #  Total TF coil peak inductive power demand, MVA
        xpower = ind_tf_total / (3600.0e0 * T_TF_CHARGE_HOURS) * itfka**2

        #  Building space:
        #  Power modules floor space, m2
        part1 = fspc1 * ntfpm * tfpmkw**0.667e0

        #  Circuit breakers floor space, m2
        part2 = fspc2 * n_tf_breakers * (v_tf_coil_dump_quench_kv * itfka) ** 0.667e0

        #  Load centres floor space, m2
        part3 = (
            fspc3 * (tfackw / (2.4e0 * nsptfc + 13.8e0 * (1.0e0 - nsptfc))) ** 0.667e0
        )

        #  Power conversion building floor area, m2
        tfcfsp = part1 + part2 + part3

        #  Dump resistor floor area, m2
        drarea = 0.5e0 * n_tf_dump_resistors * (1.0e0 + r1emj) ** 0.667e0

        #  Total TF coil power conversion building volume, m3
        tfcbv = 6.0e0 * tfcfsp

        #  TF coil AC inductive power demand, MW
        xpwrmw = xpower / 0.9e0

        #  Total steady state AC power demand, MW
        p_tf_electric_supplies_mw = (
            p_tf_electric_supplies_mw + rpower / heat_transport_variables.etatf
        )
        #  Total TF coil power conversion building floor area, m2

        # tftsp = tfcfsp
        #  Total TF coil power conversion building volume, m3

        # tftbv = tfcbv

        #  Output section
        if output:
            po.oheadr(self.outfile, "Superconducting TF Coil Power Conversion")
            po.ovarre(self.outfile, "TF coil current (kA)", "(itfka)", itfka, "OP ")
            po.ovarre(self.outfile, "Number of TF coils", "(n_tf_coils)", n_tf_coils)
            po.ovarre(
                self.outfile,
                "Voltage across a TF coil during quench (kV)",
                "(v_tf_coil_dump_quench_kv)",
                v_tf_coil_dump_quench_kv,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "TF coil charge time (hours)",
                "(T_TF_CHARGE_HOURS)",
                T_TF_CHARGE_HOURS,
            )
            po.ovarre(
                self.outfile,
                "Total inductance of TF coils (H)",
                "(ind_tf_total)",
                ind_tf_total,
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
            # po.ovarre(outfile,'Inductance per TF coil (H)','(ind_tf_coil)',ind_tf_coil, 'OP ')
            po.ovarre(self.outfile, "TF coil charging voltage (V)", "(tfcv)", tfcv)
            po.ovarre(
                self.outfile,
                "Number of DC circuit breakers",
                "(n_tf_breakers)",
                n_tf_breakers,
            )
            po.ovarre(
                self.outfile,
                "Number of dump resistors",
                "(n_tf_dump_resistors)",
                n_tf_dump_resistors,
            )
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
                self.outfile,
                "Aluminium bus current density (kA/cm2)",
                "(j_tf_bus_design_ka_cm)",
                j_tf_bus_design_ka_cm,
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
                "(m_tf_bus_aluminium_tonnes)",
                m_tf_bus_aluminium_tonnes,
                "OP ",
            )

            po.ovarre(
                self.outfile,
                "Total TF coil bus resistance (ohm)",
                "(res_tf_bus)",
                res_tf_bus,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "TF coil bus voltage drop (V)",
                "(v_tf_bus)",
                v_tf_bus,
                "OP ",
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
                "(p_tf_electric_supplies_mw)",
                p_tf_electric_supplies_mw,
                "OP ",
            )

        return (tfckw, len_tf_bus, drarea, tfcbv, p_tf_electric_supplies_mw)

    def power_profiles_over_time(
        self,
        t_precharge: float,
        t_current_ramp_up: float,
        t_fusion_ramp: float,
        t_burn: float,
        t_ramp_down: float,
        t_between_pulse: float,
        p_plant_electric_base_total_mw: float,
        p_cryo_plant_electric_mw: float,
        p_tritium_plant_electric_mw: float,
        vachtmw: float,
        p_tf_electric_supplies_mw: float,
        p_pf_electric_supplies_mw: float,
        p_coolant_pump_elec_total_mw: float,
        p_hcd_electric_total_mw: float,
        p_fusion_total_mw: float,
        p_plant_electric_gross_mw: float,
        p_plant_electric_net_mw: float,
    ) -> float:
        """
        Calculate time-dependent power profiles for different electric systems

        :param t_precharge: Precharge time (s).
        :type t_precharge: float
        :param t_current_ramp_up: Current ramp-up time (s).
        :type t_current_ramp_up: float
        :param t_fusion_ramp: Fusion ramp time (s).
        :type t_fusion_ramp: float
        :param t_burn: Burn time (s).
        :type t_burn: float
        :param t_ramp_down: Ramp-down time (s).
        :type t_ramp_down: float
        :param t_between_pulse: Time between pulses (s).
        :type t_between_pulse: float
        :param p_plant_electric_base_total_mw: Plant base electric load (MW).
        :type p_plant_electric_base_total_mw: float
        :param p_cryo_plant_electric_mw: Cryogenic plant electric load (MW).
        :type p_cryo_plant_electric_mw: float
        :param p_tritium_plant_electric_mw: Tritium plant electric load (MW).
        :type p_tritium_plant_electric_mw: float
        :param vachtmw: Vacuum pumps electric load (MW).
        :type vachtmw: float
        :param p_tf_electric_supplies_mw: TF coil electric supplies (MW).
        :type p_tf_electric_supplies_mw: float
        :param p_pf_electric_supplies_mw: PF coil electric supplies (MW).
        :type p_pf_electric_supplies_mw: float
        :param p_coolant_pump_elec_total_mw: Total coolant pump electric load (MW).
        :type p_coolant_pump_elec_total_mw: float
        :param p_hcd_electric_total_mw: HCD electric total (MW).
        :type p_hcd_electric_total_mw: float
        :param p_fusion_total_mw: Fusion power (MW).
        :type p_fusion_total_mw: float
        :param p_plant_electric_gross_mw: Gross electric power produced (MW).
        :type p_plant_electric_gross_mw: float
        :param p_plant_electric_net_mw: Net electric power produced (MW).
        :type p_plant_electric_net_mw: float

        :notes:
            - Assumes step-function changes in power at each phase transition.
            - Negative values indicate power consumption (loads).

        :returns: Total net electric energy produced over the pulse (MJ).
        :rtype: float
        """

        t_steps = np.cumsum([
            0,
            t_precharge,
            t_current_ramp_up,
            t_fusion_ramp,
            t_burn,
            t_ramp_down,
            t_between_pulse,
        ])

        # Number of time steps
        n_steps = len(t_steps)

        # Initialize arrays for each power profile
        p_fusion_total_profile_mw = np.zeros(n_steps)
        p_plant_electric_base_total_profile_mw = np.zeros(n_steps)
        p_cryo_plant_electric_profile_mw = np.zeros(n_steps)
        p_tritium_plant_electric_profile_mw = np.zeros(n_steps)
        vachtmw_profile_mw = np.zeros(n_steps)
        p_tf_electric_supplies_profile_mw = np.zeros(n_steps)
        p_pf_electric_supplies_profile_mw = np.zeros(n_steps)
        p_coolant_pump_elec_total_profile_mw = np.zeros(n_steps)
        p_hcd_electric_total_profile_mw = np.zeros(n_steps)
        p_plant_electric_gross_profile_mw = np.zeros(n_steps)
        p_plant_electric_net_profile_mw = np.zeros(n_steps)

        # Fusion power: zero until ramp-up, then during burn
        p_fusion_total_profile_mw[:2] = 0
        p_fusion_total_profile_mw[2:5] = p_fusion_total_mw
        p_fusion_total_profile_mw[5:] = 0

        # Plant base load: constant negative load throughout
        p_plant_electric_base_total_profile_mw[:] = -p_plant_electric_base_total_mw

        # Cryo plant: constant negative load throughout
        p_cryo_plant_electric_profile_mw[:] = -p_cryo_plant_electric_mw

        # Tritium plant: constant negative load throughout
        p_tritium_plant_electric_profile_mw[:] = -p_tritium_plant_electric_mw

        # Vacuum pumps: constant negative load throughout
        vachtmw_profile_mw[:] = -vachtmw

        # TF coil supplies: assume coil is always charged, so constant negative load
        p_tf_electric_supplies_profile_mw[:] = -p_tf_electric_supplies_mw

        # PF coil supplies: zero for first step, then negative during ramp-up and burn, then zero
        p_pf_electric_supplies_profile_mw[0] = 0
        p_pf_electric_supplies_profile_mw[1:5] = -p_pf_electric_supplies_mw
        p_pf_electric_supplies_profile_mw[5:] = 0

        # Coolant pump elec total: zero for first two steps, then negative during ramp-up and burn, then zero
        p_coolant_pump_elec_total_profile_mw[:2] = 0
        p_coolant_pump_elec_total_profile_mw[2:5] = -p_coolant_pump_elec_total_mw
        p_coolant_pump_elec_total_profile_mw[5:] = 0

        # HCD electric total: zero for first two steps, then negative during ramp-up and burn, then zero
        p_hcd_electric_total_profile_mw[:2] = 0
        p_hcd_electric_total_profile_mw[2:5] = -p_hcd_electric_total_mw
        p_hcd_electric_total_profile_mw[5:] = 0

        # Gross electric power: zero for first two steps, then positive during burn, then zero
        p_plant_electric_gross_profile_mw[:2] = 0
        p_plant_electric_gross_profile_mw[2:5] = p_plant_electric_gross_mw
        p_plant_electric_gross_profile_mw[5:] = 0

        # Net electric power: calculated by subtracting all loads from gross electric power
        p_plant_electric_net_profile_mw = (
            p_plant_electric_gross_profile_mw
            + p_plant_electric_base_total_profile_mw
            + p_cryo_plant_electric_profile_mw
            + p_tritium_plant_electric_profile_mw
            + vachtmw_profile_mw
            + p_tf_electric_supplies_profile_mw
            + p_pf_electric_supplies_profile_mw
            + p_coolant_pump_elec_total_profile_mw
            + p_hcd_electric_total_profile_mw
        )

        if p_plant_electric_net_profile_mw[3] != p_plant_electric_net_mw:
            logger.error(
                "Calculated net electric power during burn does not match input value."
                f"Calculated: {p_plant_electric_net_profile_mw[3]}, Input: {p_plant_electric_net_mw}"
            )

        # Integrate net electric power over the pulse to get total energy produced (MJ)
        # Assume t_steps in seconds, power in MW, so energy in MJ
        energy_made_mj = sp.integrate.trapezoid(
            p_plant_electric_net_profile_mw, t_steps
        )
        energy_made_kwh = energy_made_mj / 3.6

        return (
            energy_made_kwh,
            energy_made_mj,
            p_plant_electric_base_total_profile_mw,
            p_plant_electric_gross_profile_mw,
            p_plant_electric_net_profile_mw,
            p_hcd_electric_total_profile_mw,
            p_coolant_pump_elec_total_profile_mw,
            p_tf_electric_supplies_profile_mw,
            p_pf_electric_supplies_profile_mw,
            vachtmw_profile_mw,
            p_tritium_plant_electric_profile_mw,
            p_cryo_plant_electric_profile_mw,
            p_fusion_total_profile_mw,
        )

    def output_power_profiles_over_time(
        self,
    ):
        for i, val in enumerate(power_variables.p_plant_electric_base_total_profile_mw):
            po.ovarre(
                self.mfile,
                f"Plant total electric base load at time point {i}",
                f"(p_plant_electric_base_total_profile_mw{i})",
                val,
            )
        for i, val in enumerate(power_variables.p_plant_electric_gross_profile_mw):
            po.ovarre(
                self.mfile,
                f"Plant total electric gross at time point {i}",
                f"(p_plant_electric_gross_profile_mw{i})",
                val,
            )
        for i, val in enumerate(power_variables.p_plant_electric_net_profile_mw):
            po.ovarre(
                self.mfile,
                f"Plant total electric net at time point {i}",
                f"(p_plant_electric_net_profile_mw{i})",
                val,
            )
        for i, val in enumerate(power_variables.p_hcd_electric_total_profile_mw):
            po.ovarre(
                self.mfile,
                f"Plant total electric HCD at time point {i}",
                f"(p_hcd_electric_total_profile_mw{i})",
                val,
            )
        for i, val in enumerate(power_variables.p_coolant_pump_elec_total_profile_mw):
            po.ovarre(
                self.mfile,
                f"Plant total electric coolant pump at time point {i}",
                f"(p_coolant_pump_elec_total_profile_mw{i})",
                val,
            )
        for i, val in enumerate(power_variables.p_tf_electric_supplies_profile_mw):
            po.ovarre(
                self.mfile,
                f"Plant total electric TF supplies at time point {i}",
                f"(p_tf_electric_supplies_profile_mw{i})",
                val,
            )
        for i, val in enumerate(power_variables.p_pf_electric_supplies_profile_mw):
            po.ovarre(
                self.mfile,
                f"Plant total electric PF supplies at time point {i}",
                f"(p_pf_electric_supplies_profile_mw{i})",
                val,
            )
        for i, val in enumerate(power_variables.vachtmw_profile_mw):
            po.ovarre(
                self.mfile,
                f"Plant total electric vacuum pump power at time point {i}",
                f"(vachtmw_profile_mw{i})",
                val,
            )
        for i, val in enumerate(power_variables.p_tritium_plant_electric_profile_mw):
            po.ovarre(
                self.mfile,
                f"Plant total electric tritium plant power at time point {i}",
                f"(p_tritium_plant_electric_profile_mw{i})",
                val,
            )
        for i, val in enumerate(power_variables.p_cryo_plant_electric_profile_mw):
            po.ovarre(
                self.mfile,
                f"Plant total electric cryo plant power at time point {i}",
                f"(p_cryo_plant_electric_profile_mw{i})",
                val,
            )
        for i, val in enumerate(power_variables.p_fusion_total_profile_mw):
            po.ovarre(
                self.mfile,
                f"Plant total electric fusion plant power at time point {i}",
                f"(p_fusion_total_profile_mw{i})",
                val,
            )
