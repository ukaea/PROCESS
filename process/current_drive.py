import numpy as np

from process.fortran import (
    heat_transport_variables,
    current_drive_variables,
    physics_variables,
    cost_variables,
    constants,
    profiles_module,
    process_output as po,
    error_handling as eh,
)


class CurrentDrive:
    def __init__(self):
        self.outfile = constants.nout

    def cudriv(self, output: bool):
        """Routine to calculate the current drive power requirements
        author: P J Knight, CCFE, Culham Science Centre
        outfile : input integer : output file unit
        iprint : input integer : switch for writing to output file (1=yes)
        This routine calculates the power requirements of the current
        drive system, using a choice of models for the current drive
        efficiency.
        AEA FUS 251: A User's Guide to the PROCESS Systems Code
        """

        current_drive_variables.echpwr = 0.0e0
        current_drive_variables.pnbeam = 0.0e0
        current_drive_variables.plhybd = 0.0e0
        current_drive_variables.cnbeam = 0.0e0
        cnbeamfix = 0.0e0
        current_drive_variables.porbitlossmw = 0.0e0
        porbitlossmwfix = 0.0e0

        pinjmw1 = 0.0
        pinjmwfix = 0.0
        pinjimw1 = 0.0
        pinjemw1 = 0.0
        pinjemwfix = 0.0
        pinjimwfix = 0.0
        auxiliary_cdfix = 0.0
        faccdfix = 0.0
        gamcdfix = 0.0e0

        # To stop issues with input file we force
        # zero secondary heating if no injection method
        if current_drive_variables.iefrffix == 0:
            current_drive_variables.pheatfix = 0.0

        # check for unphysically large heating in
        # secondary injected power source
        if current_drive_variables.pheatfix > current_drive_variables.pinjfixmw:
            current_drive_variables.pheatfix = current_drive_variables.pinjfixmw

        # current_drive_variables.irfcd |  switch for current drive calculation
        # = 0   |  turned off
        # = 1   |  turned on
        if current_drive_variables.irfcd != 0:
            # put electron density in desired units (10^-20 m-3)
            dene20 = physics_variables.dene * 1.0e-20

            # If present we must calculate second current drive
            # efficiencies in units of Amps/Watt using the fixed
            # values from user input
            # current_drive_variables.iefrffix |  switch for fixed current drive efficiency model

            # Fenstermacher Lower Hybrid model
            if current_drive_variables.iefrffix == 1:
                effrfssfix = (
                    (0.36e0 * (1.0e0 + (physics_variables.te / 25.0e0) ** 1.16e0))
                    / (physics_variables.rmajor * dene20)
                    * current_drive_variables.feffcd
                )
                effcdfix = effrfssfix
            # Ion-Cyclotron current drive
            elif current_drive_variables.iefrffix == 2:
                effrfssfix = (
                    0.63e0
                    * 0.1e0
                    * physics_variables.ten
                    / (2.0e0 + physics_variables.zeff)
                    / (physics_variables.rmajor * dene20)
                    * current_drive_variables.feffcd
                )
                effcdfix = effrfssfix
            # Fenstermacher Electron Cyclotron Resonance model
            elif current_drive_variables.iefrffix == 3:
                effrfssfix = (
                    0.21e0
                    * physics_variables.ten
                    / (physics_variables.rmajor * dene20 * physics_variables.dlamee)
                    * current_drive_variables.feffcd
                )
                effcdfix = effrfssfix
            # Ehst Lower Hybrid / Fast Wave current drive
            elif current_drive_variables.iefrffix == 4:
                effrfssfix = (
                    physics_variables.te**0.77e0
                    * (0.034e0 + 0.196e0 * physics_variables.beta)
                    / (physics_variables.rmajor * dene20)
                    * (
                        32.0e0 / (5.0e0 + physics_variables.zeff)
                        + 2.0e0
                        + (12.0e0 * (6.0e0 + physics_variables.zeff))
                        / (5.0e0 + physics_variables.zeff)
                        / (3.0e0 + physics_variables.zeff)
                        + 3.76e0 / physics_variables.zeff
                    )
                    / 12.507e0
                    * current_drive_variables.feffcd
                )
                effcdfix = effrfssfix
            elif current_drive_variables.iefrffix == 5:
                (
                    effnbss,
                    current_drive_variables.fpion,
                    current_drive_variables.nbshinef,
                ) = self.iternb()
                effnbssfix = effnbss * current_drive_variables.feffcd
                effcdfix = effnbssfix
            # Culham Lower Hybrid current drive model
            elif current_drive_variables.iefrffix == 6:
                effrfss = self.cullhy()
                effrfssfix = effrfss * current_drive_variables.feffcd
                effcdfix = effrfssfix
            # Culham ECCD model
            elif current_drive_variables.iefrffix == 7:
                effrfss = self.culecd()
                effrfssfix = effrfss * current_drive_variables.feffcd
                effcdfix = effrfssfix
            # Culham Neutral Beam model
            elif current_drive_variables.iefrffix == 8:
                (
                    effnbss,
                    current_drive_variables.fpion,
                    current_drive_variables.nbshinef,
                ) = self.culnbi()
                effnbssfix = effnbss * current_drive_variables.feffcd
                effcdfix = effnbssfix
            # ECRH user input gamma
            elif current_drive_variables.iefrffix == 10:
                # Normalised current drive efficiency gamma
                current_drive_variables.gamcd = current_drive_variables.gamma_ecrh

                # Absolute current drive efficiency
                effrfssfix = current_drive_variables.gamcd / (
                    dene20 * physics_variables.rmajor
                )
                effcdfix = effrfssfix
            # EBW scaling
            elif current_drive_variables.iefrffix == 12:
                # Scaling author Simon Freethy
                # Ref : PROCESS issue 1262
                # Normalised current drive efficiency gamma
                current_drive_variables.gamcd = (
                    current_drive_variables.xi_ebw / 32.7e0
                ) * physics_variables.te

                # Absolute current drive efficiency
                effrfssfix = current_drive_variables.gamcd / (
                    dene20 * physics_variables.rmajor
                )
                effcdfix = effrfssfix

                # EBWs can only couple to plasma if cyclotron harmonic is above plasma density cut-off;
                # this behaviour is captured in the following function (ref issue #1262):
                # current_drive_variables.harnum = cyclotron harmonic number (fundamental used as default)
                # constant 'a' controls sharpness of transition
                a = 0.1e0

                fc = (
                    1.0e0
                    / (2.0e0 * np.pi)
                    * current_drive_variables.harnum
                    * constants.echarge
                    * physics_variables.bt
                    / constants.emass
                )
                fp = (
                    1.0e0
                    / (2.0e0 * np.pi)
                    * np.sqrt(
                        physics_variables.dene
                        * constants.echarge**2
                        / (constants.emass * constants.epsilon0)
                    )
                )

                density_factor = 0.5e0 * (
                    1.0e0 + np.tanh((2.0e0 / a) * ((fp - fc) / fp - a))
                )

                effcdfix = effcdfix * density_factor
                effrfssfix = effrfssfix * density_factor
            elif current_drive_variables.iefrffix == 13:
                # ECCD model for O-mode cut-off with added Te and Zeff dependance
                # Scaling author: Simon Freethy
                # Ref : PROCESS issue #2994

                fc = (
                    1
                    / (2 * np.pi)
                    * constants.echarge
                    * physics_variables.bt
                    / constants.emass
                )
                fp = (
                    1
                    / (2 * np.pi)
                    * np.sqrt(
                        (
                            (physics_variables.dene / 1.0e19)
                            * constants.echarge**2
                            / (constants.emass * constants.epsilon0)
                        )
                    )
                )

                xi_CD = 0.18e0  # This is tuned to the results of a GRAY study
                xi_CD = xi_CD * (
                    4.8e0 / (2 + physics_variables.zeff)
                )  # Zeff correction
                effrfssfix = (
                    xi_CD
                    * physics_variables.te
                    / (
                        3.27e0
                        * physics_variables.rmajor
                        * (physics_variables.dene / 1.0e19)
                    )
                )

                # O-mode case
                if current_drive_variables.wave_mode == 0:
                    f_cutoff = fp

                # X-mode case
                elif current_drive_variables.wave_mode == 1:
                    f_cutoff = 0.5 * (
                        fc
                        + np.sqrt(
                            current_drive_variables.harnum * fc**2 + 4 * fp**2
                        )
                    )

                # Plasma coupling only occurs if the plasma cut-off is below the cyclotron harmonic
                a = 0.1  # This controls how sharply the transition is reached
                cutoff_factor = 0.5 * (
                    1
                    + np.tanh(
                        (2 / (a))
                        * ((current_drive_variables.harnum * fc - f_cutoff) / fp - a)
                    )
                )
                effcdfix = effrfssfix * cutoff_factor
            elif current_drive_variables.iefrffix != 0:
                raise RuntimeError(
                    f"Current drive switch is invalid: {current_drive_variables.iefrffix = }"
                )

            if current_drive_variables.iefrffix in [1, 2, 4, 6]:
                # Injected power
                pinjemwfix = current_drive_variables.pinjfixmw

                # Wall plug power
                heat_transport_variables.pinjwpfix = (
                    current_drive_variables.pinjfixmw / current_drive_variables.etalh
                )

                # Wall plug to injector efficiency
                current_drive_variables.etacdfix = current_drive_variables.etalh

                # Normalised current drive efficiency gamma
                gamcdfix = effrfssfix * (dene20 * physics_variables.rmajor)

                # the fixed auxiliary current
                auxiliary_cdfix = (
                    effrfssfix
                    * (
                        current_drive_variables.pinjfixmw
                        - current_drive_variables.pheatfix
                    )
                    * 1.0e6
                )
                faccdfix = auxiliary_cdfix / physics_variables.plascur
            elif current_drive_variables.iefrffix in [3, 7, 10, 12, 13]:
                # Injected power
                pinjemwfix = current_drive_variables.pinjfixmw

                # Wall plug power
                heat_transport_variables.pinjwpfix = (
                    current_drive_variables.pinjfixmw / current_drive_variables.etaech
                )

                # Wall plug to injector efficiency
                current_drive_variables.etacdfix = current_drive_variables.etaech

                # the fixed auxiliary current
                auxiliary_cdfix = (
                    effrfssfix
                    * (
                        current_drive_variables.pinjfixmw
                        - current_drive_variables.pheatfix
                    )
                    * 1.0e6
                )
                faccdfix = auxiliary_cdfix / physics_variables.plascur
            elif current_drive_variables.iefrffix in [5, 8]:
                # Account for first orbit losses
                # (power due to particles that are ionised but not thermalised) [MW]:
                # This includes a second order term in shinethrough*(first orbit loss)
                current_drive_variables.forbitloss = min(
                    0.999, current_drive_variables.forbitloss
                )  # Should never be needed

                pnbitotfix = current_drive_variables.pinjfixmw / (
                    1.0e0
                    - current_drive_variables.forbitloss
                    + current_drive_variables.forbitloss
                    * current_drive_variables.nbshinef
                )

                # Shinethrough power (atoms that are not ionised) [MW]:
                nbshinemwfix = pnbitotfix * current_drive_variables.nbshinef

                # First orbit loss
                porbitlossmwfix = current_drive_variables.forbitloss * (
                    pnbitotfix - nbshinemwfix
                )

                # Power deposited
                pinjmwfix = pnbitotfix - nbshinemwfix - porbitlossmwfix
                pinjimwfix = pinjmwfix * current_drive_variables.fpion
                pinjemwfix = pinjmwfix * (1.0e0 - current_drive_variables.fpion)

                current_drive_variables.pwpnb = (
                    pnbitotfix / current_drive_variables.etanbi
                )  # neutral beam wall plug power
                heat_transport_variables.pinjwpfix = current_drive_variables.pwpnb
                current_drive_variables.etacdfix = current_drive_variables.etanbi
                gamnb = effnbssfix * (dene20 * physics_variables.rmajor)
                gamcdfix = gamnb
                cnbeamfix = (
                    1.0e-3 * (pnbitotfix * 1.0e6) / current_drive_variables.enbeam
                )  # Neutral beam current (A)
                auxiliary_cdfix = (
                    effnbssfix
                    * (
                        current_drive_variables.pinjfixmw
                        - current_drive_variables.pheatfix
                    )
                    * 1.0e6
                )
                faccdfix = auxiliary_cdfix / physics_variables.plascur

            # Fenstermacher Lower Hybrid model
            if current_drive_variables.iefrf == 1:
                effrfss = (
                    (0.36e0 * (1.0e0 + (physics_variables.te / 25.0e0) ** 1.16e0))
                    / (physics_variables.rmajor * dene20)
                    * current_drive_variables.feffcd
                )
                current_drive_variables.effcd = effrfss
            # Ion-Cyclotron current drive
            elif current_drive_variables.iefrf == 2:
                effrfss = (
                    0.63e0
                    * 0.1e0
                    * physics_variables.ten
                    / (2.0e0 + physics_variables.zeff)
                    / (physics_variables.rmajor * dene20)
                    * current_drive_variables.feffcd
                )
                current_drive_variables.effcd = effrfss
            # Fenstermacher Electron Cyclotron Resonance model
            elif current_drive_variables.iefrf == 3:
                effrfss = (
                    0.21e0
                    * physics_variables.ten
                    / (physics_variables.rmajor * dene20 * physics_variables.dlamee)
                    * current_drive_variables.feffcd
                )
                current_drive_variables.effcd = effrfss
            # Ehst Lower Hybrid / Fast Wave current drive
            elif current_drive_variables.iefrf == 4:
                effrfss = (
                    physics_variables.te**0.77e0
                    * (0.034e0 + 0.196e0 * physics_variables.beta)
                    / (physics_variables.rmajor * dene20)
                    * (
                        32.0e0 / (5.0e0 + physics_variables.zeff)
                        + 2.0e0
                        + (12.0e0 * (6.0e0 + physics_variables.zeff))
                        / (5.0e0 + physics_variables.zeff)
                        / (3.0e0 + physics_variables.zeff)
                        + 3.76e0 / physics_variables.zeff
                    )
                    / 12.507e0
                    * current_drive_variables.feffcd
                )
                current_drive_variables.effcd = effrfss
            # ITER Neutral Beam current drive
            elif current_drive_variables.iefrf == 5:
                (
                    effnbss,
                    current_drive_variables.fpion,
                    current_drive_variables.nbshinef,
                ) = self.iternb()
                effnbss = effnbss * current_drive_variables.feffcd
                current_drive_variables.effcd = effnbss
            # Culham Lower Hybrid current drive model
            elif current_drive_variables.iefrf == 6:
                effrfss = self.cullhy()
                effrfss = effrfss * current_drive_variables.feffcd
                current_drive_variables.effcd = effrfss
            # Culham ECCD model
            elif current_drive_variables.iefrf == 7:
                effrfss = self.culecd()
                effrfss = effrfss * current_drive_variables.feffcd
                current_drive_variables.effcd = effrfss
            # Culham Neutral Beam model
            elif current_drive_variables.iefrf == 8:
                (
                    effnbss,
                    current_drive_variables.fpion,
                    current_drive_variables.nbshinef,
                ) = self.culnbi()
                effnbss = effnbss * current_drive_variables.feffcd
                current_drive_variables.effcd = effnbss
            # ECRH user input gamma
            elif current_drive_variables.iefrf == 10:
                current_drive_variables.gamcd = current_drive_variables.gamma_ecrh

                # Absolute current drive efficiency
                effrfss = current_drive_variables.gamcd / (
                    dene20 * physics_variables.rmajor
                )
                current_drive_variables.effcd = effrfss
            # EBW scaling
            elif current_drive_variables.iefrf == 12:
                # Scaling author Simon Freethy
                # Ref : PROCESS issue 1262

                # Normalised current drive efficiency gamma
                current_drive_variables.gamcd = (
                    current_drive_variables.xi_ebw / 32.7e0
                ) * physics_variables.te

                # Absolute current drive efficiency
                effrfss = current_drive_variables.gamcd / (
                    dene20 * physics_variables.rmajor
                )
                current_drive_variables.effcd = effrfss
                # EBWs can only couple to plasma if cyclotron harmonic is above plasma density cut-off;
                # this behaviour is captured in the following function (ref issue #1262):
                # current_drive_variables.harnum = cyclotron harmonic number (fundamental used as default)
                # contant 'a' controls sharpness of transition
                a = 0.1e0

                fc = (
                    1.0e0
                    / (2.0e0 * np.pi)
                    * current_drive_variables.harnum
                    * constants.echarge
                    * physics_variables.bt
                    / constants.emass
                )
                fp = (
                    1.0e0
                    / (2.0e0 * np.pi)
                    * np.sqrt(
                        physics_variables.dene
                        * constants.echarge**2
                        / (constants.emass * constants.epsilon0)
                    )
                )

                density_factor = 0.5e0 * (
                    1.0e0 + np.tanh((2.0e0 / a) * ((fp - fc) / fp - a))
                )

                current_drive_variables.effcd = (
                    current_drive_variables.effcd * density_factor
                )
                effrfss = effrfss * density_factor

            elif current_drive_variables.iefrf == 13:
                # ECCD model for O-mode cut-off with added Te and Zeff dependance
                # Scaling author: Simon Freethy
                # Ref : PROCESS issue #2994

                fc = (
                    1
                    / (2 * np.pi)
                    * constants.echarge
                    * physics_variables.bt
                    / constants.emass
                )
                fp = (
                    1
                    / (2 * np.pi)
                    * np.sqrt(
                        (
                            (physics_variables.dene / 1.0e19)
                            * constants.echarge**2
                            / (constants.emass * constants.epsilon0)
                        )
                    )
                )

                xi_CD = 0.18e0  # This is tuned to the results of a GRAY study
                xi_CD = xi_CD * (
                    4.8e0 / (2 + physics_variables.zeff)
                )  # Zeff correction
                effrfss = (
                    xi_CD
                    * physics_variables.te
                    / (
                        3.27e0
                        * physics_variables.rmajor
                        * (physics_variables.dene / 1.0e19)
                    )
                )

                # O-mode case
                if current_drive_variables.wave_mode == 0:
                    f_cutoff = fp

                # X-mode case
                elif current_drive_variables.wave_mode == 1:
                    f_cutoff = 0.5 * (
                        fc
                        + np.sqrt(
                            current_drive_variables.harnum * fc**2 + 4 * fp**2
                        )
                    )

                # Plasma coupling only occurs if the plasma cut-off is below the cyclotron harmonic
                a = 0.1  # This controls how sharply the transition is reached
                cutoff_factor = 0.5 * (
                    1
                    + np.tanh(
                        (2 / (a))
                        * ((current_drive_variables.harnum * fc - f_cutoff) / fp - a)
                    )
                )
                current_drive_variables.effcd = effrfss * cutoff_factor
            else:
                raise RuntimeError(
                    f"Current drive switch is invalid: {current_drive_variables.iefrf = }"
                )

            # Compute current drive wall plug and injected powers (MW) and efficiencies
            auxiliary_cd = physics_variables.faccd * physics_variables.plascur

            # LHCD or ICCD
            if current_drive_variables.iefrf in [1, 2, 4, 6]:
                # Injected power
                current_drive_variables.plhybd = (
                    1.0e-6
                    * (physics_variables.faccd - faccdfix)
                    * physics_variables.plascur
                    / effrfss
                    + current_drive_variables.pheat
                )
                pinjimw1 = 0.0e0
                pinjemw1 = current_drive_variables.plhybd

                # Wall plug power
                current_drive_variables.pwplh = (
                    current_drive_variables.plhybd / current_drive_variables.etalh
                )
                pinjwp1 = current_drive_variables.pwplh

                # Wall plug to injector efficiency
                current_drive_variables.etacd = current_drive_variables.etalh

                # Normalised current drive efficiency gamma
                gamrf = effrfss * (dene20 * physics_variables.rmajor)
                current_drive_variables.gamcd = gamrf
            # ECCD
            elif current_drive_variables.iefrf in [3, 7, 10, 12, 13]:
                # Injected power (set to close to close the Steady-state current equilibrium)
                current_drive_variables.echpwr = (
                    1.0e-6
                    * (physics_variables.faccd - faccdfix)
                    * physics_variables.plascur
                    / effrfss
                    + current_drive_variables.pheat
                )
                pinjemw1 = current_drive_variables.echpwr

                # Wall plug power
                current_drive_variables.echwpow = (
                    current_drive_variables.echpwr / current_drive_variables.etaech
                )

                # Wall plug to injector efficiency
                pinjwp1 = current_drive_variables.echwpow
                current_drive_variables.etacd = current_drive_variables.etaech
            elif current_drive_variables.iefrf in [5, 8]:
                # MDK. See Gitlab issue #248, and scanned note.
                power1 = (
                    1.0e-6
                    * (physics_variables.faccd - faccdfix)
                    * physics_variables.plascur
                    / effnbss
                    + current_drive_variables.pheat
                )

                # Account for first orbit losses
                # (power due to particles that are ionised but not thermalised) [MW]:
                # This includes a second order term in shinethrough*(first orbit loss)
                current_drive_variables.forbitloss = min(
                    0.999, current_drive_variables.forbitloss
                )  # Should never be needed

                current_drive_variables.pnbitot = power1 / (
                    1.0e0
                    - current_drive_variables.forbitloss
                    + current_drive_variables.forbitloss
                    * current_drive_variables.nbshinef
                )

                # Shinethrough power (atoms that are not ionised) [MW]:
                current_drive_variables.nbshinemw = (
                    current_drive_variables.pnbitot * current_drive_variables.nbshinef
                )

                # First orbit loss
                current_drive_variables.porbitlossmw = (
                    current_drive_variables.forbitloss
                    * (
                        current_drive_variables.pnbitot
                        - current_drive_variables.nbshinemw
                    )
                )

                # Power deposited
                pinjmw1 = (
                    current_drive_variables.pnbitot
                    - current_drive_variables.nbshinemw
                    - current_drive_variables.porbitlossmw
                )
                pinjimw1 = pinjmw1 * current_drive_variables.fpion
                pinjemw1 = pinjmw1 * (1.0e0 - current_drive_variables.fpion)

                current_drive_variables.pwpnb = (
                    current_drive_variables.pnbitot / current_drive_variables.etanbi
                )  # neutral beam wall plug power
                pinjwp1 = current_drive_variables.pwpnb
                current_drive_variables.etacd = current_drive_variables.etanbi
                gamnb = effnbss * (dene20 * physics_variables.rmajor)
                current_drive_variables.gamcd = gamnb
                current_drive_variables.cnbeam = (
                    1.0e-3
                    * (current_drive_variables.pnbitot * 1.0e6)
                    / current_drive_variables.enbeam
                )  # Neutral beam current (A)

            # Total injected power
            # sum contributions from primary and secondary systems
            current_drive_variables.pinjmw = (
                pinjemw1 + pinjimw1 + pinjemwfix + pinjimwfix
            )
            pinjmw1 = pinjemw1 + pinjimw1
            pinjmwfix = pinjemwfix + pinjimwfix
            current_drive_variables.pinjemw = pinjemw1 + pinjemwfix
            current_drive_variables.pinjimw = pinjimw1 + pinjimwfix
            heat_transport_variables.pinjwp = (
                pinjwp1 + heat_transport_variables.pinjwpfix
            )

            # Reset injected power to zero for ignited plasma (fudge)
            if physics_variables.ignite == 1:
                heat_transport_variables.pinjwp = 0.0e0

            # Ratio of fusion to input (injection+ohmic) power
            if (
                abs(
                    current_drive_variables.pinjmw
                    + current_drive_variables.porbitlossmw
                    + physics_variables.pohmmw
                )
                < 1.0e-6
            ):
                current_drive_variables.bigq = 1.0e18
            else:
                current_drive_variables.bigq = physics_variables.powfmw / (
                    current_drive_variables.pinjmw
                    + current_drive_variables.porbitlossmw
                    + physics_variables.pohmmw
                )

        if not output:
            return

        po.oheadr(self.outfile, "Current Drive System")

        if current_drive_variables.irfcd == 0:
            po.ocmmnt(self.outfile, "No current drive used")
            po.oblnkl(self.outfile)
            return

        if current_drive_variables.iefrf in [1, 4, 6]:
            po.ocmmnt(self.outfile, "Lower Hybrid Current Drive")
        elif current_drive_variables.iefrf == 2:
            po.ocmmnt(self.outfile, "Ion Cyclotron Current Drive")
        elif current_drive_variables.iefrf in [3, 7]:
            po.ocmmnt(self.outfile, "Electron Cyclotron Current Drive")
        elif current_drive_variables.iefrf in [5, 8]:
            po.ocmmnt(self.outfile, "Neutral Beam Current Drive")
        elif current_drive_variables.iefrf == 10:
            po.ocmmnt(
                self.outfile, "Electron Cyclotron Current Drive (user input gamma_CD)"
            )
        elif current_drive_variables.iefrf == 12:
            po.ocmmnt(self.outfile, "EBW current drive")
        elif current_drive_variables.iefrf == 13:
            po.ocmmnt(
                self.outfile,
                "Electron Cyclotron Current Drive (O-mode cutoff with Zeff & Te)",
            )

        po.ovarin(
            self.outfile,
            "Current drive efficiency model",
            "(iefrf)",
            current_drive_variables.iefrf,
        )

        if current_drive_variables.iefrffix in [1, 4, 6]:
            po.ocmmnt(self.outfile, "Lower Hybrid Current Drive")
        elif current_drive_variables.iefrffix == 2:
            po.ocmmnt(self.outfile, "Ion Cyclotron Current Drive")
        elif current_drive_variables.iefrffix in [3, 7]:
            po.ocmmnt(self.outfile, "Electron Cyclotron Current Drive")
        elif current_drive_variables.iefrffix in [5, 8]:
            po.ocmmnt(self.outfile, "Neutral Beam Current Drive")
        elif current_drive_variables.iefrffix == 10:
            po.ocmmnt(
                self.outfile, "Electron Cyclotron Current Drive (user input gamma_CD)"
            )
        elif current_drive_variables.iefrffix == 12:
            po.ocmmnt(self.outfile, "EBW current drive")
        elif current_drive_variables.iefrffix == 13:
            po.ocmmnt(
                self.outfile,
                "Electron Cyclotron Current Drive (O-mode cutoff with Zeff & Te)",
            )

        po.ovarin(
            self.outfile,
            "Secondary current drive efficiency model",
            "(iefrffix)",
            current_drive_variables.iefrffix,
        )

        if physics_variables.ignite == 1:
            po.ocmmnt(
                self.outfile,
                "Ignited plasma; injected power only used for start-up phase",
            )

        po.oblnkl(self.outfile)

        if abs(physics_variables.facoh) > 1.0e-8:
            po.ocmmnt(self.outfile, "Current is driven by both inductive")
            po.ocmmnt(self.outfile, "and non-inductive means.")

        po.ovarre(
            self.outfile,
            "Ratio of power for flat-top to start-up (MW)",
            "(startupratio)",
            cost_variables.startupratio,
        )
        po.ovarre(
            self.outfile,
            "Auxiliary power used for plasma heating only (MW)",
            "(pheat)",
            current_drive_variables.pheat + current_drive_variables.pheatfix,
        )
        po.ovarre(
            self.outfile,
            "Power injected for current drive (MW)",
            "(pcurrentdrivemw)",
            current_drive_variables.pinjmw
            - current_drive_variables.pheat
            - current_drive_variables.pheatfix,
        )
        po.ovarre(
            self.outfile,
            "Maximum Allowed Bootstrap current fraction",
            "(bscfmax)",
            current_drive_variables.bscfmax,
        )
        if current_drive_variables.iefrffix != 0:
            po.ovarre(
                self.outfile,
                "Power injected for main current drive (MW)",
                "(pcurrentdrivemw1)",
                pinjmw1 - current_drive_variables.pheat,
            )
            po.ovarre(
                self.outfile,
                "Power injected for secondary current drive (MW)",
                "(pcurrentdrivemw2)",
                pinjmwfix - current_drive_variables.pheatfix,
            )

        po.ovarre(
            self.outfile,
            "Fusion gain factor Q",
            "(bigq)",
            current_drive_variables.bigq,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Auxiliary current drive (A)",
            "(auxiliary_cd)",
            auxiliary_cd,
            "OP ",
        )
        if current_drive_variables.iefrffix != 0:
            po.ovarre(
                self.outfile,
                "Secondary auxiliary current drive (A)",
                "(auxiliary_cdfix)",
                auxiliary_cdfix,
                "OP ",
            )

        po.ovarre(
            self.outfile,
            "Current drive efficiency (A/W)",
            "(effcd)",
            current_drive_variables.effcd,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Normalised current drive efficiency, gamma (10^20 A/W-m2)",
            "(gamcd)",
            current_drive_variables.gamcd,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Wall plug to injector efficiency",
            "(etacd)",
            current_drive_variables.etacd,
        )

        if current_drive_variables.iefrf == 10:
            po.ovarre(
                self.outfile,
                "ECRH plasma heating efficiency",
                "(gamma_ecrh)",
                current_drive_variables.gamma_ecrh,
            )
        if current_drive_variables.iefrf == 12:
            po.ovarre(
                self.outfile,
                "EBW plasma heating efficiency",
                "(xi_ebw)",
                current_drive_variables.xi_ebw,
            )
        if current_drive_variables.iefrf in [12, 13]:
            po.ovarre(
                self.outfile,
                "EC harmonic number",
                "(harnum)",
                current_drive_variables.harnum,
            )
        if current_drive_variables.iefrf == 13:
            po.ovarin(
                self.outfile,
                "EC cutoff wave mode switch",
                "(wave_mode)",
                current_drive_variables.wave_mode,
            )

        if current_drive_variables.iefrffix != 0:
            po.ovarre(
                self.outfile,
                "Secondary current drive efficiency (A/W)",
                "(effcdfix)",
                effcdfix,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Seconday wall plug to injector efficiency",
                "(etacdfix)",
                current_drive_variables.etacdfix,
            )
            po.ovarre(
                self.outfile,
                "Normalised secondary current drive efficiency, gamma (10^20 A/W-m2)",
                "(gamcdfix)",
                gamcdfix,
                "OP ",
            )

        po.osubhd(self.outfile, "Fractions of current drive :")
        po.ovarrf(
            self.outfile,
            "Bootstrap fraction",
            "(bootipf)",
            current_drive_variables.bootipf,
            "OP ",
        )
        po.ovarrf(
            self.outfile,
            "Diamagnetic fraction",
            "(diaipf)",
            current_drive_variables.diaipf,
            "OP ",
        )
        po.ovarrf(
            self.outfile,
            "Pfirsch-Schlueter fraction",
            "(psipf)",
            current_drive_variables.psipf,
            "OP ",
        )
        po.ovarrf(
            self.outfile,
            "Auxiliary current drive fraction",
            "(faccd)",
            physics_variables.faccd,
            "OP ",
        )
        po.ovarrf(
            self.outfile,
            "Inductive fraction",
            "(facoh)",
            physics_variables.facoh,
            "OP ",
        )
        # Add total error check.
        po.ovarrf(
            self.outfile,
            "Total",
            "(plasipf+faccd+facoh)",
            current_drive_variables.plasipf
            + physics_variables.faccd
            + physics_variables.facoh,
        )
        if (
            abs(
                current_drive_variables.plasipf
                + physics_variables.faccd
                + physics_variables.facoh
                - 1.0e0
            )
            > 1.0e-8
        ):
            po.ocmmnt(self.outfile, "ERROR: current drive fractions do not add to 1")

        # MDK Add physics_variables.fvsbrnni as it can be an iteration variable
        po.ovarrf(
            self.outfile,
            "Fraction of the plasma current produced by non-inductive means",
            "(fvsbrnni)",
            physics_variables.fvsbrnni,
        )

        if (
            abs(current_drive_variables.bootipf - current_drive_variables.bscfmax)
            < 1.0e-8
        ):
            po.ocmmnt(self.outfile, "Warning : bootstrap current fraction is at")
            po.ocmmnt(self.outfile, "          its prescribed maximum.")

        po.oblnkl(self.outfile)

        if abs(current_drive_variables.plhybd) > 1.0e-8:
            po.ovarre(self.outfile, "RF efficiency (A/W)", "(effrfss)", effrfss, "OP ")
            po.ovarre(self.outfile, "RF gamma (10^20 A/W-m2)", "(gamrf)", gamrf, "OP ")
            po.ovarre(
                self.outfile,
                "Lower hybrid injected power (MW)",
                "(plhybd)",
                current_drive_variables.plhybd,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Lower hybrid wall plug efficiency",
                "(etalh)",
                current_drive_variables.etalh,
            )
            po.ovarre(
                self.outfile,
                "Lower hybrid wall plug power (MW)",
                "(pwplh)",
                current_drive_variables.pwplh,
                "OP ",
            )

        # MDK rearranged and added current_drive_variables.nbshinemw
        # if (abs(current_drive_variables.pnbeam) > 1.0e-8) :
        if (
            (current_drive_variables.iefrf == 5)
            or (current_drive_variables.iefrf == 8)
            or (current_drive_variables.iefrffix == 5)
            or (current_drive_variables.iefrffix == 8)
        ):
            po.ovarre(
                self.outfile,
                "Neutral beam energy (keV)",
                "(enbeam)",
                current_drive_variables.enbeam,
            )
            if (current_drive_variables.iefrf == 5) or (
                current_drive_variables.iefrf == 8
            ):
                po.ovarre(
                    self.outfile,
                    "Neutral beam current (A)",
                    "(cnbeam)",
                    current_drive_variables.cnbeam,
                    "OP ",
                )

            if (current_drive_variables.iefrffix == 5) or (
                current_drive_variables.iefrffix == 8
            ):
                po.ovarre(
                    self.outfile,
                    "Secondary fixed neutral beam current (A)",
                    "(cnbeamfix)",
                    cnbeamfix,
                    "OP ",
                )

            if (current_drive_variables.iefrf == 5) or (
                current_drive_variables.iefrf == 8
            ):
                po.ovarre(
                    self.outfile, "Beam efficiency (A/W)", "(effnbss)", effnbss, "OP "
                )

            if (current_drive_variables.iefrffix == 5) or (
                current_drive_variables.iefrffix == 8
            ):
                po.ovarre(
                    self.outfile,
                    "Secondary fixed beam efficiency (A/W)",
                    "(effnbssfix)",
                    effnbssfix,
                    "OP ",
                )

            po.ovarre(
                self.outfile, "Beam gamma (10^20 A/W-m2)", "(gamnb)", gamnb, "OP "
            )
            po.ovarre(
                self.outfile,
                "Neutral beam wall plug efficiency",
                "(etanbi)",
                current_drive_variables.etanbi,
            )
            po.ovarre(
                self.outfile,
                "Beam decay lengths to centre",
                "(taubeam)",
                current_drive_variables.taubeam,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Beam shine-through fraction",
                "(nbshinef)",
                current_drive_variables.nbshinef,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Neutral beam wall plug power (MW)",
                "(pwpnb)",
                current_drive_variables.pwpnb,
                "OP ",
            )

            po.oblnkl(self.outfile)
            po.ocmmnt(self.outfile, "Neutral beam power balance :")
            po.ocmmnt(self.outfile, "----------------------------")
            if (current_drive_variables.iefrf == 5) or (
                current_drive_variables.iefrf == 8
            ):
                po.ovarrf(
                    self.outfile,
                    "Beam first orbit loss power (MW)",
                    "(porbitlossmw)",
                    current_drive_variables.porbitlossmw,
                    "OP ",
                )
                po.ovarrf(
                    self.outfile,
                    "Beam shine-through power [MW]",
                    "(nbshinemw)",
                    current_drive_variables.nbshinemw,
                    "OP ",
                )
                po.ovarrf(
                    self.outfile,
                    "Beam power deposited in plasma (MW)",
                    "(pinjmw)",
                    pinjmw1,
                    "OP ",
                )
                po.ovarrf(
                    self.outfile,
                    "Maximum allowable beam power (MW)",
                    "(pinjalw)",
                    current_drive_variables.pinjalw,
                )
                po.ovarrf(
                    self.outfile,
                    "Total (MW)",
                    "(current_drive_variables.porbitlossmw+current_drive_variables.nbshinemw+current_drive_variables.pinjmw)",
                    current_drive_variables.porbitlossmw
                    + current_drive_variables.nbshinemw
                    + pinjmw1,
                )
                po.oblnkl(self.outfile)
                po.ovarrf(
                    self.outfile,
                    "Beam power entering vacuum vessel (MW)",
                    "(pnbitot)",
                    current_drive_variables.pnbitot,
                    "OP ",
                )

            if (current_drive_variables.iefrffix == 5) or (
                current_drive_variables.iefrffix == 8
            ):
                po.oblnkl(self.outfile)
                po.ocmmnt(self.outfile, "Secondary fixed neutral beam power balance :")
                po.ocmmnt(self.outfile, "----------------------------")
                po.ovarrf(
                    self.outfile,
                    "Secondary fixed beam first orbit loss power (MW)",
                    "(porbitlossmwfix)",
                    porbitlossmwfix,
                    "OP ",
                )
                po.ovarrf(
                    self.outfile,
                    "Secondary fixed beam shine-through power [MW]",
                    "(nbshinemwfix)",
                    nbshinemwfix,
                    "OP ",
                )
                po.ovarrf(
                    self.outfile,
                    "Secondary fixed beam power deposited in plasma (MW)",
                    "(pinjmwfix)",
                    pinjmwfix,
                    "OP ",
                )
                po.ovarrf(
                    self.outfile,
                    "Maximum allowable beam power (MW)",
                    "(pinjalw)",
                    current_drive_variables.pinjalw,
                )
                po.ovarrf(
                    self.outfile,
                    "Secondary fixed total (MW)",
                    "(porbitlossmwfixed+nbshinemwfix+pinjmwfix)",
                    porbitlossmwfix + nbshinemwfix + pinjmwfix,
                )
                po.oblnkl(self.outfile)
                po.ovarrf(
                    self.outfile,
                    "Secondary beam power entering vacuum vessel (MW)",
                    "(pnbitotfix)",
                    pnbitotfix,
                    "OP ",
                )

            po.oblnkl(self.outfile)

            po.ovarre(
                self.outfile,
                "Fraction of beam energy to ions",
                "(fpion)",
                current_drive_variables.fpion,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Beam duct shielding thickness (m)",
                "(nbshield)",
                current_drive_variables.nbshield,
            )
            po.ovarre(
                self.outfile,
                "Beam tangency radius / Plasma major radius",
                "(frbeam)",
                current_drive_variables.frbeam,
            )
            po.ovarre(
                self.outfile,
                "Beam centreline tangency radius (m)",
                "(rtanbeam)",
                current_drive_variables.rtanbeam,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Maximum possible tangency radius (m)",
                "(rtanmax)",
                current_drive_variables.rtanmax,
                "OP ",
            )

        if abs(current_drive_variables.echpwr) > 1.0e-8:
            po.ovarre(
                self.outfile,
                "Electron cyclotron injected power (MW)",
                "(echpwr)",
                current_drive_variables.echpwr,
                "OP ",
            )
            po.ovarrf(
                self.outfile,
                "Maximum allowable ECRH power (MW)",
                "(pinjalw)",
                current_drive_variables.pinjalw,
            )
            po.ovarre(
                self.outfile,
                "ECH wall plug efficiency",
                "(etaech)",
                current_drive_variables.etaech,
            )
            po.ovarre(
                self.outfile,
                "ECH wall plug power (MW)",
                "(echwpow)",
                current_drive_variables.echwpow,
                "OP ",
            )

        if abs(current_drive_variables.pinjfixmw) > 1.0e-8:
            po.ovarrf(
                self.outfile,
                "Fixed ECRH power (MW)",
                "(pinjmwfix)",
                current_drive_variables.pinjmwfix,
            )
            po.ovarre(
                self.outfile,
                "ECH wall plug efficiency",
                "(etaech)",
                current_drive_variables.etaech,
            )
            po.ovarre(
                self.outfile,
                "Secondary fixed ECH wall plug power (MW)",
                "(pinjwpfix)",
                current_drive_variables.pinjwpfix,
                "OP ",
            )

    def iternb(self):
        """Routine to calculate ITER Neutral Beam current drive parameters
        author: P J Knight, CCFE, Culham Science Centre
        effnbss : output real : neutral beam current drive efficiency (A/W)
        fpion   : output real : fraction of NB power given to ions
        fshine  : output real : shine-through fraction of beam
        This routine calculates the current drive parameters for a
        neutral beam system, based on the 1990 ITER model.
        AEA FUS 251: A User's Guide to the PROCESS Systems Code
        ITER Physics Design Guidelines: 1989 [IPDG89], N. A. Uckan et al,
        ITER Documentation Series No.10, IAEA/ITER/DS/10, IAEA, Vienna, 1990
        """
        # Check argument sanity
        if (1 + physics_variables.eps) < current_drive_variables.frbeam:
            eh.fdiags[0] = physics_variables.eps
            eh.fdiags[1] = current_drive_variables.frbeam
            eh.report_error(15)

        # Calculate beam path length to centre
        dpath = physics_variables.rmajor * np.sqrt(
            (1.0 + physics_variables.eps) ** 2 - current_drive_variables.frbeam**2
        )

        # Calculate beam stopping cross-section
        sigstop = self.sigbeam(
            current_drive_variables.enbeam / physics_variables.abeam,
            physics_variables.te,
            physics_variables.dene,
            physics_variables.ralpne,
            physics_variables.rncne,
            physics_variables.rnone,
            physics_variables.rnfene,
        )

        # Calculate number of decay lengths to centre
        current_drive_variables.taubeam = dpath * physics_variables.dene * sigstop

        # Shine-through fraction of beam
        fshine = np.exp(-2.0 * dpath * physics_variables.dene * sigstop)
        fshine = max(fshine, 1e-20)

        # Deuterium and tritium beam densities
        dend = physics_variables.deni * (1.0 - current_drive_variables.ftritbm)
        dent = physics_variables.deni * current_drive_variables.ftritbm

        # Power split to ions / electrons
        fpion = self.cfnbi(
            physics_variables.abeam,
            current_drive_variables.enbeam,
            physics_variables.ten,
            physics_variables.dene,
            dend,
            dent,
            physics_variables.zeffai,
            physics_variables.dlamie,
        )

        # Current drive efficiency
        effnbss = current_drive_variables.frbeam * self.etanb(
            physics_variables.abeam,
            physics_variables.alphan,
            physics_variables.alphat,
            physics_variables.aspect,
            physics_variables.dene,
            current_drive_variables.enbeam,
            physics_variables.rmajor,
            physics_variables.ten,
            physics_variables.zeff,
        )

        return effnbss, fpion, fshine

    def cullhy(self):
        """Routine to calculate Lower Hybrid current drive efficiency
        author: P J Knight, CCFE, Culham Science Centre
        effrfss : output real : lower hybrid current drive efficiency (A/W)
        This routine calculates the current drive parameters for a
        lower hybrid system, based on the AEA FUS 172 model.
        AEA FUS 251: A User's Guide to the PROCESS Systems Code
        AEA FUS 172: Physics Assessment for the European Reactor Study
        """
        rratio = self.lhrad()
        rpenet = rratio * physics_variables.rminor

        # Local density, temperature, toroidal field at this minor radius

        dlocal = 1.0e-19 * profiles_module.nprofile(
            rratio,
            physics_variables.rhopedn,
            physics_variables.ne0,
            physics_variables.neped,
            physics_variables.nesep,
            physics_variables.alphan,
        )
        tlocal = profiles_module.tprofile(
            rratio,
            physics_variables.rhopedt,
            physics_variables.te0,
            physics_variables.teped,
            physics_variables.tesep,
            physics_variables.alphat,
            physics_variables.tbeta,
        )
        blocal = (
            physics_variables.bt
            * physics_variables.rmajor
            / (physics_variables.rmajor - rpenet)
        )  # Calculated on inboard side

        # Parallel refractive index needed for plasma access

        frac = np.sqrt(dlocal) / blocal
        nplacc = frac + np.sqrt(1.0e0 + frac * frac)

        # Local inverse aspect ratio

        epslh = rpenet / physics_variables.rmajor

        # LH normalised efficiency (A/W m**-2)

        x = 24.0e0 / (nplacc * np.sqrt(tlocal))

        term01 = 6.1e0 / (nplacc * nplacc * (physics_variables.zeff + 5.0e0))
        term02 = 1.0e0 + (tlocal / 25.0e0) ** 1.16e0
        term03 = epslh**0.77e0 * np.sqrt(12.25e0 + x * x)
        term04 = 3.5e0 * epslh**0.77e0 + x

        if term03 > term04:
            eh.fdiags[0] = term03
            eh.fdiags[1] = term04
            eh.report_error(129)

        gamlh = term01 * term02 * (1.0e0 - term03 / term04)

        # Current drive efficiency (A/W)

        return gamlh / ((0.1e0 * dlocal) * physics_variables.rmajor)

    def culecd(self):
        """Routine to calculate Electron Cyclotron current drive efficiency
        author: M R O'Brien, CCFE, Culham Science Centre
        author: P J Knight, CCFE, Culham Science Centre
        effrfss : output real : electron cyclotron current drive efficiency (A/W)
        This routine calculates the current drive parameters for a
        electron cyclotron system, based on the AEA FUS 172 model.
        AEA FUS 251: A User's Guide to the PROCESS Systems Code
        AEA FUS 172: Physics Assessment for the European Reactor Study
        """
        rrr = 1.0e0 / 3.0e0

        #  Temperature
        tlocal = profiles_module.tprofile(
            rrr,
            physics_variables.rhopedt,
            physics_variables.te0,
            physics_variables.teped,
            physics_variables.tesep,
            physics_variables.alphat,
            physics_variables.tbeta,
        )

        #  Density (10**20 m**-3)
        dlocal = 1.0e-20 * profiles_module.nprofile(
            rrr,
            physics_variables.rhopedn,
            physics_variables.ne0,
            physics_variables.neped,
            physics_variables.nesep,
            physics_variables.alphan,
        )

        #  Inverse aspect ratio
        epsloc = rrr * physics_variables.rminor / physics_variables.rmajor

        #  Effective charge (use average value)
        zlocal = physics_variables.zeff

        #  Coulomb logarithm for ion-electron collisions
        #  (From J. A. Wesson, 'Tokamaks', Clarendon Press, Oxford, p.293)
        coulog = 15.2e0 - 0.5e0 * np.log(dlocal) + np.log(tlocal)

        #  Calculate normalised current drive efficiency at four different
        #  poloidal angles, and average.
        #  cosang = cosine of the poloidal angle at which ECCD takes place
        #         = +1 outside, -1 inside.
        cosang = 1.0e0
        ecgam1 = self.eccdef(tlocal, epsloc, zlocal, cosang, coulog)
        cosang = 0.5e0
        ecgam2 = self.eccdef(tlocal, epsloc, zlocal, cosang, coulog)
        cosang = -0.5e0
        ecgam3 = self.eccdef(tlocal, epsloc, zlocal, cosang, coulog)
        cosang = -1.0e0
        ecgam4 = self.eccdef(tlocal, epsloc, zlocal, cosang, coulog)

        #  Normalised current drive efficiency (A/W m**-2)
        ecgam = 0.25e0 * (ecgam1 + ecgam2 + ecgam3 + ecgam4)

        #  Current drive efficiency (A/W)
        return ecgam / (dlocal * physics_variables.rmajor)

    def eccdef(self, tlocal, epsloc, zlocal, cosang, coulog):
        """Routine to calculate Electron Cyclotron current drive efficiency
        author: M R O'Brien, CCFE, Culham Science Centre
        author: P J Knight, CCFE, Culham Science Centre
        tlocal : input real : local electron temperature (keV)
        epsloc : input real : local inverse aspect ratio
        zlocal : input real : local plasma effective charge
        cosang : input real : cosine of the poloidal angle at which ECCD takes
        place (+1 outside, -1 inside)
        coulog : input real : local coulomb logarithm for ion-electron collisions
        ecgam  : output real : normalised current drive efficiency (A/W m**-2)
        This routine calculates the current drive parameters for a
        electron cyclotron system, based on the AEA FUS 172 model.
        It works out the ECCD efficiency using the formula due to Cohen
        quoted in the ITER Physics Design Guidelines : 1989
        (but including division by the Coulomb Logarithm omitted from
        IPDG89). We have assumed gamma**2-1 << 1, where gamma is the
        relativistic factor. The notation follows that in IPDG89.
        <P>The answer ECGAM is the normalised efficiency nIR/P with n the
        local density in 10**20 /m**3, I the driven current in MAmps,
        R the major radius in metres, and P the absorbed power in MWatts.
        AEA FUS 251: A User's Guide to the PROCESS Systems Code
        AEA FUS 172: Physics Assessment for the European Reactor Study
        ITER Physics Design Guidelines: 1989 [IPDG89], N. A. Uckan et al,
        ITER Documentation Series No.10, IAEA/ITER/DS/10, IAEA, Vienna, 1990
        """
        mcsq = 9.1095e-31 * 2.9979e8**2 / (1.0e3 * 1.6022e-19)  # keV
        f = 16.0e0 * (tlocal / mcsq) ** 2

        #  fp is the derivative of f with respect to gamma, the relativistic
        #  factor, taken equal to 1 + 2T/(m c**2)

        fp = 16.0e0 * tlocal / mcsq

        #  lam is IPDG89's lambda. LEGEND calculates the Legendre function of
        #  order alpha and argument lam, palpha, and its derivative, palphap.
        #  Here alpha satisfies alpha(alpha+1) = -8/(1+zlocal). alpha is of the
        #  form  (-1/2 + ix), with x a real number and i = sqrt(-1).

        lam = 1.0e0
        palpha, palphap = self.legend(zlocal, lam)

        lams = np.sqrt(2.0e0 * epsloc / (1.0e0 + epsloc))
        palphas, _ = self.legend(zlocal, lams)

        #  hp is the derivative of IPDG89's h function with respect to lam

        h = -4.0e0 * lam / (zlocal + 5.0e0) * (1.0e0 - lams * palpha / (lam * palphas))
        hp = -4.0e0 / (zlocal + 5.0e0) * (1.0e0 - lams * palphap / palphas)

        #  facm is IPDG89's momentum conserving factor

        facm = 1.5e0
        y = mcsq / (2.0e0 * tlocal) * (1.0e0 + epsloc * cosang)

        #  We take the negative of the IPDG89 expression to get a positive
        #  number

        ecgam = (
            -7.8e0
            * facm
            * np.sqrt((1.0e0 + epsloc) / (1.0e0 - epsloc))
            / coulog
            * (h * fp - 0.5e0 * y * f * hp)
        )

        if ecgam < 0.0e0:
            eh.report_error(17)
        return ecgam

    def culnbi(self):
        """Routine to calculate Neutral Beam current drive parameters
        author: P J Knight, CCFE, Culham Science Centre
        effnbss : output real : neutral beam current drive efficiency (A/W)
        fpion   : output real : fraction of NB power given to ions
        fshine  : output real : shine-through fraction of beam
        This routine calculates Neutral Beam current drive parameters
        using the corrections outlined in AEA FUS 172 to the ITER method.
        <P>The result cannot be guaranteed for devices with aspect ratios far
        from that of ITER (approx. 2.8).
        AEA FUS 251: A User's Guide to the PROCESS Systems Code
        AEA FUS 172: Physics Assessment for the European Reactor Study
        """
        if (1.0e0 + physics_variables.eps) < current_drive_variables.frbeam:
            eh.fdiags[0] = physics_variables.eps
            eh.fdiags[1] = current_drive_variables.frbeam
            eh.report_error(20)

        #  Calculate beam path length to centre

        dpath = physics_variables.rmajor * np.sqrt(
            (1.0e0 + physics_variables.eps) ** 2 - current_drive_variables.frbeam**2
        )

        #  Calculate beam stopping cross-section

        sigstop = self.sigbeam(
            current_drive_variables.enbeam / physics_variables.abeam,
            physics_variables.te,
            physics_variables.dene,
            physics_variables.ralpne,
            physics_variables.rncne,
            physics_variables.rnone,
            physics_variables.rnfene,
        )

        #  Calculate number of decay lengths to centre

        current_drive_variables.taubeam = dpath * physics_variables.dnla * sigstop

        #  Shine-through fraction of beam

        fshine = np.exp(-2.0e0 * dpath * physics_variables.dnla * sigstop)
        fshine = max(fshine, 1.0e-20)

        #  Deuterium and tritium beam densities

        dend = physics_variables.deni * (1.0e0 - current_drive_variables.ftritbm)
        dent = physics_variables.deni * current_drive_variables.ftritbm

        #  Power split to ions / electrons

        fpion = self.cfnbi(
            physics_variables.abeam,
            current_drive_variables.enbeam,
            physics_variables.ten,
            physics_variables.dene,
            dend,
            dent,
            physics_variables.zeffai,
            physics_variables.dlamie,
        )

        #  Current drive efficiency

        effnbss = self.etanb2(
            physics_variables.abeam,
            physics_variables.alphan,
            physics_variables.alphat,
            physics_variables.aspect,
            physics_variables.dene,
            physics_variables.dnla,
            current_drive_variables.enbeam,
            current_drive_variables.frbeam,
            fshine,
            physics_variables.rmajor,
            physics_variables.rminor,
            physics_variables.ten,
            physics_variables.zeff,
        )

        return effnbss, fpion, fshine

    def lhrad(self):
        """Routine to calculate Lower Hybrid wave absorption radius
        author: P J Knight, CCFE, Culham Science Centre
        rratio  : output real : minor radius of penetration / rminor
        This routine determines numerically the minor radius at which the
        damping of Lower Hybrid waves occurs, using a Newton-Raphson method.
        AEA FUS 251: A User's Guide to the PROCESS Systems Code
        AEA FUS 172: Physics Assessment for the European Reactor Study
        """
        #  Correction to refractive index (kept within valid bounds)
        drfind = min(0.7e0, max(0.1e0, 12.5e0 / physics_variables.te0))

        #  Use Newton-Raphson method to establish the correct minor radius
        #  ratio. g is calculated as a function of r / r_minor, where g is
        #  the difference between the results of the two formulae for the
        #  energy E given in AEA FUS 172, p.58. The required minor radius
        #  ratio has been found when g is sufficiently close to zero.

        #  Initial guess for the minor radius ratio

        rat0 = 0.8e0

        for _ in range(100):
            #  Minor radius ratios either side of the latest guess

            r1 = rat0 - 1.0e-3 * rat0
            r2 = rat0 + 1.0e-3 * rat0

            #  Evaluate g at rat0, r1, r2

            g0 = self.lheval(drfind, rat0)
            g1 = self.lheval(drfind, r1)
            g2 = self.lheval(drfind, r2)

            #  Calculate gradient of g with respect to minor radius ratio

            dgdr = (g2 - g1) / (r2 - r1)

            #  New approximation

            rat1 = rat0 - g0 / dgdr

            #  Force this approximation to lie within bounds

            rat1 = max(0.0001e0, rat1)
            rat1 = min(0.9999e0, rat1)

            if abs(g0) <= 0.01e0:
                break
            rat0 = rat1

        else:
            eh.report_error(16)
            rat0 = 0.8e0

        return rat0

    def etanb2(
        self,
        abeam,
        alphan,
        alphat,
        aspect,
        dene,
        dnla,
        enbeam,
        frbeam,
        fshine,
        rmajor,
        rminor,
        ten,
        zeff,
    ):
        """Routine to find neutral beam current drive efficiency
        using the ITER 1990 formulation, plus correction terms
        outlined in Culham Report AEA FUS 172
        author: P J Knight, CCFE, Culham Science Centre
        abeam   : input real : beam ion mass (amu)
        alphan  : input real : density profile factor
        alphat  : input real : temperature profile factor
        aspect  : input real : aspect ratio
        dene    : input real : volume averaged electron density (m**-3)
        dnla    : input real : line averaged electron density (m**-3)
        enbeam  : input real : neutral beam energy (keV)
        frbeam  : input real : R_tangent / R_major for neutral beam injection
        fshine  : input real : shine-through fraction of beam
        rmajor  : input real : plasma major radius (m)
        rminor  : input real : plasma minor radius (m)
        ten     : input real : density weighted average electron temperature (keV)
        zeff    : input real : plasma effective charge
        This routine calculates the current drive efficiency in A/W of
        a neutral beam system, based on the 1990 ITER model,
        plus correction terms outlined in Culham Report AEA FUS 172.
        <P>The formulae are from AEA FUS 172, unless denoted by IPDG89.
        AEA FUS 251: A User's Guide to the PROCESS Systems Code
        AEA FUS 172: Physics Assessment for the European Reactor Study
        ITER Physics Design Guidelines: 1989 [IPDG89], N. A. Uckan et al,
        ITER Documentation Series No.10, IAEA/ITER/DS/10, IAEA, Vienna, 1990
        """
        #  Charge of beam ions
        zbeam = 1.0

        #  Fitting factor (IPDG89)
        bbd = 1.0

        #  Volume averaged electron density (10**20 m**-3)
        dene20 = dene / 1e20

        #  Line averaged electron density (10**20 m**-3)
        dnla20 = dnla / 1e20

        #  Critical energy (MeV) (power to electrons = power to ions) (IPDG89)
        #  N.B. ten is in keV
        ecrit = 0.01 * abeam * ten

        #  Beam energy in MeV
        ebmev = enbeam / 1e3

        #  x and y coefficients of function J0(x,y) (IPDG89)
        xjs = ebmev / (bbd * ecrit)
        xj = np.sqrt(xjs)

        yj = 0.8 * zeff / abeam

        #  Fitting function J0(x,y)
        j0 = xjs / (4.0 + 3.0 * yj + xjs * (xj + 1.39 + 0.61 * yj**0.7))

        #  Effective inverse aspect ratio, with a limit on its maximum value
        epseff = min(0.2, (0.5 / aspect))

        #  Reduction in the reverse electron current
        #  due to neoclassical effects
        gfac = (1.55 + 0.85 / zeff) * np.sqrt(epseff) - (0.2 + 1.55 / zeff) * epseff

        # Reduction in the net beam driven current
        #  due to the reverse electron current
        ffac = 1.0 - (zbeam / zeff) * (1.0 - gfac)

        #  Normalisation to allow results to be valid for
        #  non-ITER plasma size and density:

        #  Line averaged electron density (10**20 m**-3) normalised to ITER
        nnorm = 1.0

        #  Distance along beam to plasma centre
        r = max(rmajor, rmajor * frbeam)
        eps1 = rminor / r

        if (1.0 + eps1) < frbeam:
            eh.fdiags[0] = eps1
            eh.fdiags[1] = frbeam
            eh.report_error(21)

        d = rmajor * np.sqrt((1.0 + eps1) ** 2 - frbeam**2)

        # Distance along beam to plasma centre for ITER
        # assuming a tangency radius equal to the major radius
        epsitr = 2.15 / 6.0
        dnorm = 6.0 * np.sqrt(2.0 * epsitr + epsitr**2)

        #  Normalisation to beam energy (assumes a simplified formula for
        #  the beam stopping cross-section)
        ebnorm = ebmev * ((nnorm * dnorm) / (dnla20 * d)) ** (1.0 / 0.78)

        #  A_bd fitting coefficient, after normalisation with ebnorm
        abd = (
            0.107
            * (1.0 - 0.35 * alphan + 0.14 * alphan**2)
            * (1.0 - 0.21 * alphat)
            * (1.0 - 0.2 * ebnorm + 0.09 * ebnorm**2)
        )

        #  Normalised current drive efficiency (A/W m**-2) (IPDG89)
        gamnb = 5.0 * abd * 0.1 * ten * (1.0 - fshine) * frbeam * j0 / 0.2 * ffac

        #  Current drive efficiency (A/W)
        return gamnb / (dene20 * rmajor)

    def lheval(self, drfind, rratio):
        """Routine to evaluate the difference between electron energy
        expressions required to find the Lower Hybrid absorption radius
        author: P J Knight, CCFE, Culham Science Centre
        drfind  : input real : correction to parallel refractive index
        rratio  : input real : guess for radius of penetration / rminor
        ediff   : output real : difference between the E values (keV)
        This routine evaluates the difference between the values calculated
        from the two equations for the electron energy E, given in
        AEA FUS 172, p.58. This difference is used to locate the Lower Hybrid
        wave absorption radius via a Newton-Raphson method, in calling
        routine <A HREF="lhrad.html">lhrad</A>.
        AEA FUS 251: A User's Guide to the PROCESS Systems Code
        AEA FUS 172: Physics Assessment for the European Reactor Study
        """
        dlocal = 1.0e-19 * profiles_module.nprofile(
            rratio,
            physics_variables.rhopedn,
            physics_variables.ne0,
            physics_variables.neped,
            physics_variables.nesep,
            physics_variables.alphan,
        )

        #  Local electron temperature

        tlocal = profiles_module.tprofile(
            rratio,
            physics_variables.rhopedt,
            physics_variables.te0,
            physics_variables.teped,
            physics_variables.tesep,
            physics_variables.alphat,
            physics_variables.tbeta,
        )

        #  Local toroidal field (evaluated at the inboard region of the flux surface)

        blocal = (
            physics_variables.bt
            * physics_variables.rmajor
            / (physics_variables.rmajor - rratio * physics_variables.rminor)
        )

        #  Parallel refractive index needed for plasma access

        frac = np.sqrt(dlocal) / blocal
        nplacc = frac + np.sqrt(1.0e0 + frac * frac)

        #  Total parallel refractive index

        refind = nplacc + drfind

        #  First equation for electron energy E

        e1 = 511.0e0 * (np.sqrt(1.0e0 + 1.0e0 / (refind * refind)) - 1.0e0)

        #  Second equation for E

        e2 = 7.0e0 * tlocal

        #  Difference

        return e1 - e2

    def etanb(self, abeam, alphan, alphat, aspect, dene, ebeam, rmajor, ten, zeff):
        """Routine to find neutral beam current drive efficiency
        using the ITER 1990 formulation
        author: P J Knight, CCFE, Culham Science Centre
        abeam   : input real : beam ion mass (amu)
        alphan  : input real : density profile factor
        alphat  : input real : temperature profile factor
        aspect  : input real : aspect ratio
        dene    : input real : volume averaged electron density (m**-3)
        enbeam  : input real : neutral beam energy (keV)
        rmajor  : input real : plasma major radius (m)
        ten     : input real : density weighted average electron temp. (keV)
        zeff    : input real : plasma effective charge
        This routine calculates the current drive efficiency of
        a neutral beam system, based on the 1990 ITER model.
        AEA FUS 251: A User's Guide to the PROCESS Systems Code
        ITER Physics Design Guidelines: 1989 [IPDG89], N. A. Uckan et al,
        ITER Documentation Series No.10, IAEA/ITER/DS/10, IAEA, Vienna, 1990
        """

        zbeam = 1.0
        bbd = 1.0

        dene20 = 1e-20 * dene

        # Ratio of E_beam/E_crit
        xjs = ebeam / (bbd * 10.0 * abeam * ten)
        xj = np.sqrt(xjs)

        yj = 0.8 * zeff / abeam

        rjfunc = xjs / (4.0 + 3.0 * yj + xjs * (xj + 1.39 + 0.61 * yj**0.7))

        epseff = 0.5 / aspect
        gfac = (1.55 + 0.85 / zeff) * np.sqrt(epseff) - (0.2 + 1.55 / zeff) * epseff
        ffac = 1.0 / zbeam - (1.0 - gfac) / zeff

        abd = (
            0.107
            * (1.0 - 0.35 * alphan + 0.14 * alphan**2)
            * (1.0 - 0.21 * alphat)
            * (1.0 - 0.2e-3 * ebeam + 0.09e-6 * ebeam**2)
        )

        return abd * (5.0 / rmajor) * (0.1 * ten / dene20) * rjfunc / 0.2 * ffac

    def cfnbi(self, afast, efast, te, ne, nd, nt, zeffai, xlmbda):
        """Routine to calculate the fraction of the fast particle energy
        coupled to the ions
        author: P J Knight, CCFE, Culham Science Centre
        afast   : input real : mass of fast particle (units of proton mass)
        efast   : input real : energy of fast particle (keV)
        te      : input real : density weighted average electron temp. (keV)
        ne      : input real : volume averaged electron density (m**-3)
        nd      : input real : deuterium beam density (m**-3)
        nt      : input real : tritium beam density (m**-3)
        zeffai  : input real : mass weighted plasma effective charge
        xlmbda  : input real : ion-electron coulomb logarithm
        fpion   : output real : fraction of fast particle energy coupled to ions
        This routine calculates the fast particle energy coupled to
        the ions in the neutral beam system.
        AEA FUS 251: A User's Guide to the PROCESS Systems Code
        """
        # atmd = 2.0
        atmdt = 2.5
        # atmt = 3.0
        c = 3.0e8
        me = 9.1e-31
        # zd = 1.0
        # zt = 1.0

        # xlbd = self.xlmbdabi(afast, atmd, efast, te, ne)
        # xlbt = self.xlmbdabi(afast, atmt, efast, te, ne)

        # sum = nd * zd * zd * xlbd / atmd + nt * zt * zt * xlbt / atmt
        # ecritfix = 16.0e0 * te * afast * (sum / (ne * xlmbda)) ** (2.0e0 / 3.0e0)

        xlmbdai = self.xlmbdabi(afast, atmdt, efast, te, ne)
        sumln = zeffai * xlmbdai / xlmbda
        xlnrat = (3.0e0 * np.sqrt(np.pi) / 4.0e0 * me / constants.mproton * sumln) ** (
            2.0e0 / 3.0e0
        )
        ve = c * np.sqrt(2.0e0 * te / 511.0e0)

        ecritfi = (
            afast
            * constants.mproton
            * ve
            * ve
            * xlnrat
            / (2.0e0 * constants.echarge * 1.0e3)
        )

        x = np.sqrt(efast / ecritfi)
        t1 = np.log((x * x - x + 1.0e0) / ((x + 1.0e0) ** 2))
        thx = (2.0e0 * x - 1.0e0) / np.sqrt(3.0e0)
        t2 = 2.0e0 * np.sqrt(3.0e0) * (np.atan(thx) + np.pi / 6.0e0)

        return (t1 + t2) / (3.0e0 * x * x)

    def xlmbdabi(self, mb, mth, eb, t, nelec):
        """Calculates the Coulomb logarithm for ion-ion collisions
        author: P J Knight, CCFE, Culham Science Centre
        mb     : input real : mass of fast particle (units of proton mass)
        mth    : input real : mass of background ions (units of proton mass)
        eb     : input real : energy of fast particle (keV)
        t      : input real : density weighted average electron temp. (keV)
        nelec  : input real : volume averaged electron density (m**-3)
        This function calculates the Coulomb logarithm for ion-ion
        collisions where the relative velocity may be large compared
        with the background ('mt') thermal velocity.
        Mikkelson and Singer, Nuc Tech/Fus, 4, 237 (1983)
        """

        x1 = (t / 10.0) * (eb / 1000.0) * mb / (nelec / 1e20)
        x2 = mth / (mth + mb)

        return 23.7 + np.log(x2 * np.sqrt(x1))

    def legend(self, zlocal, arg):
        """Routine to calculate Legendre function and its derivative
        author: M R O'Brien, CCFE, Culham Science Centre
        author: P J Knight, CCFE, Culham Science Centre
        zlocal  : input real : local plasma effective charge
        arg     : input real : argument of Legendre function
        palpha  : output real : value of Legendre function
        palphap : output real : derivative of Legendre function
        This routine calculates the Legendre function <CODE>palpha</CODE>
        of argument <CODE>arg</CODE> and order
        <CODE>alpha = -0.5 + i sqrt(xisq)</CODE>,
        and its derivative <CODE>palphap</CODE>.
        <P>This Legendre function is a conical function and we use the series
        in <CODE>xisq</CODE> given in Abramowitz and Stegun. The
        derivative is calculated from the derivative of this series.
        <P>The derivatives were checked by calculating <CODE>palpha</CODE> for
        neighbouring arguments. The calculation of <CODE>palpha</CODE> for zero
        argument was checked by comparison with the expression
        <CODE>palpha(0) = 1/sqrt(pi) * cos(pi*alpha/2) * gam1 / gam2</CODE>
        (Abramowitz and Stegun, eqn 8.6.1). Here <CODE>gam1</CODE> and
        <CODE>gam2</CODE> are the Gamma functions of arguments
        <CODE>0.5*(1+alpha)</CODE> and <CODE>0.5*(2+alpha)</CODE> respectively.
        Abramowitz and Stegun, equation 8.12.1
        """
        if abs(arg) > (1.0e0 + 1.0e-10):
            eh.fdiags[0] = arg
            eh.report_error(18)

        arg2 = min(arg, (1.0e0 - 1.0e-10))
        sinsq = 0.5e0 * (1.0e0 - arg2)
        xisq = 0.25e0 * (32.0e0 * zlocal / (zlocal + 1.0e0) - 1.0e0)
        palpha = 1.0e0
        pold = 1.0e0
        pterm = 1.0e0
        palphap = 0.0e0
        poldp = 0.0e0

        for n in range(10000):
            #  Check for convergence every 20 iterations

            if (n > 1) and ((n % 20) == 1):
                term1 = 1.0e-10 * max(abs(pold), abs(palpha))
                term2 = 1.0e-10 * max(abs(poldp), abs(palphap))

                if (abs(pold - palpha) < term1) and (abs(poldp - palphap) < term2):
                    return palpha, palphap

                pold = palpha
                poldp = palphap

            pterm = (
                pterm
                * (4.0e0 * xisq + (2.0e0 * n - 1.0e0) ** 2)
                / (2.0e0 * n) ** 2
                * sinsq
            )
            palpha = palpha + pterm
            palphap = palphap - n * pterm / (1.0e0 - arg2)
        else:
            eh.report_error(19)

    def sigbeam(self, eb, te, ne, rnhe, rnc, rno, rnfe):
        """Calculates the stopping cross-section for a hydrogen
        beam in a fusion plasma
        author: P J Knight, CCFE, Culham Science Centre
        eb     : input real : beam energy (kev/amu)
        te     : input real : electron temperature (keV)
        ne     : input real : electron density (10^20m-3)
        rnhe   : input real : alpha density / ne
        rnc    : input real : carbon density /ne
        rno    : input real : oxygen density /ne
        rnfe   : input real : iron density /ne
        This function calculates the stopping cross-section (m^-2)
        for a hydrogen beam in a fusion plasma.
        Janev, Boley and Post, Nuclear Fusion 29 (1989) 2125
        """
        a = np.array(
            [
                [
                    [4.4, -2.49e-2],
                    [7.46e-2, 2.27e-3],
                    [3.16e-3, -2.78e-5],
                ],
                [
                    [2.3e-1, -1.15e-2],
                    [-2.55e-3, -6.2e-4],
                    [1.32e-3, 3.38e-5],
                ],
            ]
        )

        b = np.array(
            [
                [
                    [[-2.36, -1.49, -1.41, -1.03], [0.185, -0.0154, -4.08e-4, 0.106]],
                    [
                        [-0.25, -0.119, -0.108, -0.0558],
                        [-0.0381, -0.015, -0.0138, -3.72e-3],
                    ],
                ],
                [
                    [
                        [0.849, 0.518, 0.477, 0.322],
                        [-0.0478, 7.18e-3, 1.57e-3, -0.0375],
                    ],
                    [
                        [0.0677, 0.0292, 0.0259, 0.0124],
                        [0.0105, 3.66e-3, 3.33e-3, 8.61e-4],
                    ],
                ],
                [
                    [
                        [-0.0588, -0.0336, -0.0305, -0.0187],
                        [4.34e-3, 3.41e-4, 7.35e-4, 3.53e-3],
                    ],
                    [
                        [-4.48e-3, -1.79e-3, -1.57e-3, -7.43e-4],
                        [-6.76e-4, -2.04e-4, -1.86e-4, -5.12e-5],
                    ],
                ],
            ]
        )

        z = np.array([2.0, 6.0, 8.0, 26.0])
        nn = np.array([rnhe, rnc, rno, rnfe])

        nen = ne * 1e-19

        s1 = 0.0
        for k in range(2):
            for j in range(3):
                for i in range(2):
                    s1 += (
                        a[i, j, k]
                        * (np.log(eb)) ** i
                        * (np.log(nen)) ** j
                        * (np.log(te)) ** k
                    )

        sz = 0.0
        for l in range(4):  # noqa: E741
            for k in range(2):
                for j in range(2):
                    for i in range(3):
                        sz += (
                            b[i, j, k, l]
                            * (np.log(eb)) ** i
                            * (np.log(nen)) ** j
                            * (np.log(te)) ** k
                            * nn[l]
                            * z[l]
                            * (z[l] - 1.0)
                        )

        return max(1e-20 * (np.exp(s1) / eb * (1.0 + sz)), 1e-23)
