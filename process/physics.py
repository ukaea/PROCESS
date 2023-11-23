import numpy
import scipy.integrate as integrate
import math
import process.physics_functions as physics_funcs
from process.utilities.f2py_string_patch import f2py_compatible_to_string
from process.fortran import (
    current_drive_module,
    constraint_variables,
    reinke_variables,
    reinke_module,
    impurity_radiation_module,
    constants,
    physics_functions_module,
    physics_variables,
    physics_module,
    pulse_variables,
    times_variables,
    current_drive_variables,
    error_handling,
    fwbs_variables,
    build_variables,
    divertor_variables,
    numerics,
    stellarator_variables,
    process_output as po,
    profiles_module,
)


class Physics:
    def __init__(self, plasma_profile):
        self.outfile = constants.nout
        self.plasma_profile = plasma_profile

    def physics(self):
        """
        Routine to calculate tokamak plasma physics information
        author: P J Knight, CCFE, Culham Science Centre
        None
        This routine calculates all the primary plasma physics
        M. Kovari et al, 2014, "PROCESS": A systems code for fusion power plants -
        Part 1: Physics https://www.sciencedirect.com/science/article/pii/S0920379614005961
        H. Zohm et al, 2013, On the Physics Guidelines for a Tokamak DEMO
        https://iopscience.iop.org/article/10.1088/0029-5515/53/7/073019
        T. Hartmann, 2013, Development of a modular systems code to analyse the
        implications of physics assumptions on the design of a demonstration fusion power plant
        https://inis.iaea.org/search/search.aspx?orig_q=RN:45031642
        """
        # kappaa_IPB = physics_variables.vol / (
        #     2.0e0
        #     * numpy.pi
        #     * numpy.pi
        #     * physics_variables.rminor
        #     * physics_variables.rminor
        #     * physics_variables.rmajor
        # )

        if physics_variables.icurr == 2:
            physics_variables.q95 = (
                physics_variables.q * 1.3e0 * (1.0e0 - physics_variables.eps) ** 0.6e0
            )
        else:
            physics_variables.q95 = (
                physics_variables.q
            )  # i.e. input (or iteration variable) value

        # Calculate plasma composition
        # Issue #261 Remove old radiation model (imprad_model=0)
        physics_module.plasma_composition()

        #  Calculate plasma current
        (
            physics_variables.bp,
            physics_variables.qstar,
            physics_variables.plascur,
        ) = physics_module.culcur(
            physics_variables.alphaj,
            physics_variables.alphap,
            physics_variables.bt,
            physics_variables.eps,
            physics_variables.icurr,
            physics_variables.iprofile,
            physics_variables.kappa,
            physics_variables.kappa95,
            physics_variables.p0,
            physics_variables.pperim,
            physics_variables.q0,
            physics_variables.q,
            physics_variables.rli,
            physics_variables.rmajor,
            physics_variables.rminor,
            physics_variables.sf,
            physics_variables.triang,
            physics_variables.triang95,
        )

        #  Calculate density and temperature profile quantities
        #  If physics_variables.ipedestal = 1 then set pedestal density to
        #    physics_variables.fgwped * Greenwald density limit
        #  Note: this used to be done before plasma current
        if (physics_variables.ipedestal == 1) and (physics_variables.fgwped >= 0e0):
            physics_variables.neped = (
                physics_variables.fgwped
                * 1.0e14
                * physics_variables.plascur
                / (numpy.pi * physics_variables.rminor * physics_variables.rminor)
            )

        if (physics_variables.ipedestal == 1) and (physics_variables.fgwsep >= 0e0):
            physics_variables.nesep = (
                physics_variables.fgwsep
                * 1.0e14
                * physics_variables.plascur
                / (numpy.pi * physics_variables.rminor * physics_variables.rminor)
            )

        self.plasma_profile.run()

        # Calculate total magnetic field [T]
        physics_variables.btot = physics_functions_module.total_mag_field()

        # Calculate physics_variables.beta poloidal [-]
        physics_variables.betap = physics_functions_module.beta_poloidal()

        #  Set PF coil ramp times
        if pulse_variables.lpulse != 1:
            if times_variables.tohsin == 0.0e0:
                times_variables.tohs = physics_variables.plascur / 5.0e5
                times_variables.tramp = times_variables.tohs
                times_variables.tqnch = times_variables.tohs
            else:
                times_variables.tohs = times_variables.tohsin

        else:
            if times_variables.pulsetimings == 0.0e0:
                # times_variables.tramp is input
                times_variables.tohs = physics_variables.plascur / 1.0e5
                times_variables.tqnch = times_variables.tohs

            else:
                #  times_variables.tohs is set either in INITIAL or INPUT, or by being
                #  iterated using limit equation 41.
                times_variables.tramp = max(times_variables.tramp, times_variables.tohs)
                # tqnch = max(tqnch,tohs)
                times_variables.tqnch = times_variables.tohs

        #  Reset second times_variables.tburn value (times_variables.tburn0).
        #  This is used to ensure that the burn time is used consistently;
        #  see convergence loop in fcnvmc1, evaluators.f90
        times_variables.tburn0 = times_variables.tburn

        #  Pulse and down times : The reactor is assumed to be 'down'
        #  at all times outside of the plasma current flat-top period.
        #  The pulse length is the duration of non-zero plasma current
        times_variables.tpulse = (
            times_variables.tohs
            + times_variables.theat
            + times_variables.tburn
            + times_variables.tqnch
        )
        times_variables.tdown = (
            times_variables.tramp
            + times_variables.tohs
            + times_variables.tqnch
            + times_variables.tdwell
        )

        #  Total cycle time
        times_variables.tcycle = (
            times_variables.tramp
            + times_variables.tohs
            + times_variables.theat
            + times_variables.tburn
            + times_variables.tqnch
            + times_variables.tdwell
        )

        #  Calculate bootstrap current fraction using various models
        current_drive_variables.bscf_iter89 = self.bootstrap_fraction_iter89(
            physics_variables.aspect,
            physics_variables.beta,
            physics_variables.btot,
            current_drive_variables.cboot,
            physics_variables.plascur,
            physics_variables.q95,
            physics_variables.q0,
            physics_variables.rmajor,
            physics_variables.vol,
        )

        # Profile parameters are meaningless with ipedestal=3
        betat = (
            physics_variables.beta
            * physics_variables.btot**2
            / physics_variables.bt**2
        )
        current_drive_variables.bscf_nevins = (
            current_drive_variables.cboot
            * self.bootstrap_fraction_nevins(
                physics_variables.alphan,
                physics_variables.alphat,
                betat,
                physics_variables.bt,
                physics_variables.dene,
                physics_variables.plascur,
                physics_variables.q95,
                physics_variables.q0,
                physics_variables.rmajor,
                physics_variables.rminor,
                physics_variables.ten,
                physics_variables.zeff,
            )
        )

        #  Wilson scaling uses thermal poloidal beta, not total
        betpth = (
            physics_variables.beta - physics_variables.betaft - physics_variables.betanb
        ) * (physics_variables.btot / physics_variables.bp) ** 2
        current_drive_variables.bscf_wilson = (
            current_drive_variables.cboot
            * self.bootstrap_fraction_wilson(
                physics_variables.alphaj,
                physics_variables.alphap,
                physics_variables.alphat,
                betpth,
                physics_variables.q0,
                physics_variables.q95,
                physics_variables.rmajor,
                physics_variables.rminor,
            )
        )

        # Hender scaling for diamagnetic current at tight physics_variables.aspect ratio
        current_drive_variables.diacf_hender = self.diamagnetic_fraction_hender(
            physics_variables.beta
        )

        # SCENE scaling for diamagnetic current
        current_drive_variables.diacf_scene = self.diamagnetic_fraction_scene(
            physics_variables.beta, physics_variables.q95, physics_variables.q0
        )

        # Pfirsch-Schlüter scaling for diamagnetic current
        current_drive_variables.pscf_scene = self.ps_fraction_scene(
            physics_variables.beta
        )

        current_drive_variables.bscf_sauter = (
            current_drive_variables.cboot * self.bootstrap_fraction_sauter()
        )

        if current_drive_variables.bscfmax < 0.0e0:
            current_drive_variables.bootipf = abs(current_drive_variables.bscfmax)
            current_drive_variables.plasipf = current_drive_variables.bootipf
        else:
            if physics_variables.ibss == 1:
                current_drive_variables.bootipf = current_drive_variables.bscf_iter89
            elif physics_variables.ibss == 2:
                current_drive_variables.bootipf = current_drive_variables.bscf_nevins
            elif physics_variables.ibss == 3:
                current_drive_variables.bootipf = current_drive_variables.bscf_wilson
            elif physics_variables.ibss == 4:
                current_drive_variables.bootipf = current_drive_variables.bscf_sauter
            else:
                error_handling.idiags[0] = physics_variables.ibss
                error_handling.report_error(75)

            physics_module.err242 = 0
            if current_drive_variables.bootipf > current_drive_variables.bscfmax:
                current_drive_variables.bootipf = min(
                    current_drive_variables.bootipf, current_drive_variables.bscfmax
                )
                physics_module.err242 = 1

            if physics_variables.idia == 1:
                current_drive_variables.diaipf = current_drive_variables.diacf_hender
            elif physics_variables.idia == 2:
                current_drive_variables.diaipf = current_drive_variables.diacf_scene

            if physics_variables.ips == 1:
                current_drive_variables.psipf = current_drive_variables.pscf_scene

            current_drive_variables.plasipf = (
                current_drive_variables.bootipf
                + current_drive_variables.diaipf
                + current_drive_variables.psipf
            )

        #  Plasma driven current fraction (Bootstrap + Diamagnetic
        #  + Pfirsch-Schlüter) constrained to be less than
        #  or equal to the total fraction of the plasma current
        #  produced by non-inductive means (which also includes
        #  the current drive proportion)
        physics_module.err243 = 0
        if current_drive_variables.plasipf > physics_variables.fvsbrnni:
            current_drive_variables.plasipf = min(
                current_drive_variables.plasipf, physics_variables.fvsbrnni
            )
            physics_module.err243 = 1

        #  Fraction of plasma current produced by inductive means
        physics_variables.facoh = max(1.0e-10, (1.0e0 - physics_variables.fvsbrnni))
        #   Fraction of plasma current produced by auxiliary current drive
        physics_variables.faccd = (
            physics_variables.fvsbrnni - current_drive_variables.plasipf
        )

        #  Auxiliary current drive power calculations

        if current_drive_variables.irfcd != 0:
            current_drive_module.cudriv(constants.nout, 0)

        # Calculate fusion power + components

        # physics_funcs.palph(self.plasma_profile)
        fusion_rate = physics_funcs.FusionReactionRate(self.plasma_profile)
        fusion_rate.calculate_fusion_rates()
        fusion_rate.set_physics_variables()

        #
        physics_variables.pdt = physics_module.pdtpv * physics_variables.vol
        physics_variables.pdhe3 = physics_module.pdhe3pv * physics_variables.vol
        physics_variables.pdd = physics_module.pddpv * physics_variables.vol

        #  Calculate neutral beam slowing down effects
        #  If ignited, then ignore beam fusion effects

        if (current_drive_variables.cnbeam != 0.0e0) and (
            physics_variables.ignite == 0
        ):
            (
                physics_variables.betanb,
                physics_variables.dnbeam2,
                physics_variables.palpnb,
            ) = physics_functions_module.beamfus(
                physics_variables.beamfus0,
                physics_variables.betbm0,
                physics_variables.bp,
                physics_variables.bt,
                current_drive_variables.cnbeam,
                physics_variables.dene,
                physics_variables.deni,
                physics_variables.dlamie,
                physics_variables.ealphadt,
                current_drive_variables.enbeam,
                physics_variables.fdeut,
                physics_variables.ftrit,
                current_drive_variables.ftritbm,
                physics_module.sigvdt,
                physics_variables.ten,
                physics_variables.tin,
                physics_variables.vol,
                physics_variables.zeffai,
            )
            physics_variables.fusionrate = (
                physics_variables.fusionrate
                + 1.0e6
                * physics_variables.palpnb
                / (1.0e3 * physics_variables.ealphadt * constants.echarge)
                / physics_variables.vol
            )
            physics_variables.alpharate = (
                physics_variables.alpharate
                + 1.0e6
                * physics_variables.palpnb
                / (1.0e3 * physics_variables.ealphadt * constants.echarge)
                / physics_variables.vol
            )

        physics_variables.pdt = physics_variables.pdt + 5.0e0 * physics_variables.palpnb

        # Create some derived values and add beam contribution to fusion power
        (
            physics_variables.palpmw,
            physics_variables.pneutmw,
            physics_variables.pchargemw,
            physics_variables.betaft,
            physics_variables.palpipv,
            physics_variables.palpepv,
            physics_variables.pfuscmw,
            physics_variables.powfmw,
        ) = physics_functions_module.palph2(
            physics_variables.bt,
            physics_variables.bp,
            physics_variables.dene,
            physics_variables.deni,
            physics_variables.dnitot,
            physics_variables.falpe,
            physics_variables.falpi,
            physics_variables.palpnb,
            physics_variables.ifalphap,
            physics_variables.pchargepv,
            physics_variables.pneutpv,
            physics_variables.ten,
            physics_variables.tin,
            physics_variables.vol,
            physics_variables.palppv,
        )

        #  Nominal mean neutron wall load on entire first wall area including divertor and beam holes
        #  Note that 'fwarea' excludes these, so they have been added back in.
        if physics_variables.iwalld == 1:
            physics_variables.wallmw = (
                physics_variables.ffwal
                * physics_variables.pneutmw
                / physics_variables.sarea
            )
        else:
            if physics_variables.idivrt == 2:
                # Double null configuration
                physics_variables.wallmw = (
                    (1.0e0 - fwbs_variables.fhcd - 2.0e0 * fwbs_variables.fdiv)
                    * physics_variables.pneutmw
                    / build_variables.fwarea
                )
            else:
                # Single null Configuration
                physics_variables.wallmw = (
                    (1.0e0 - fwbs_variables.fhcd - fwbs_variables.fdiv)
                    * physics_variables.pneutmw
                    / build_variables.fwarea
                )

        #  Calculate ion/electron equilibration power

        physics_variables.piepv = physics_module.rether(
            physics_variables.alphan,
            physics_variables.alphat,
            physics_variables.dene,
            physics_variables.dlamie,
            physics_variables.te,
            physics_variables.ti,
            physics_variables.zeffai,
        )

        #  Calculate radiation power

        radpwrdata = physics_funcs.radpwr(self.plasma_profile)
        physics_variables.pbrempv = radpwrdata.pbrempv
        physics_variables.plinepv = radpwrdata.plinepv
        physics_variables.psyncpv = radpwrdata.psyncpv
        physics_variables.pcoreradpv = radpwrdata.pcoreradpv
        physics_variables.pedgeradpv = radpwrdata.pedgeradpv
        physics_variables.pradpv = radpwrdata.pradpv

        physics_variables.pinnerzoneradmw = (
            physics_variables.pcoreradpv * physics_variables.vol
        )
        physics_variables.pouterzoneradmw = (
            physics_variables.pedgeradpv * physics_variables.vol
        )
        physics_variables.pradmw = physics_variables.pradpv * physics_variables.vol

        #  Calculate ohmic power
        (
            physics_variables.pohmpv,
            physics_variables.pohmmw,
            physics_variables.rpfac,
            physics_variables.rplas,
        ) = physics_module.pohm(
            physics_variables.facoh,
            physics_variables.kappa95,
            physics_variables.plascur,
            physics_variables.rmajor,
            physics_variables.rminor,
            physics_variables.ten,
            physics_variables.vol,
            physics_variables.zeff,
        )

        #  Calculate L- to H-mode power threshold for different scalings

        physics_variables.pthrmw = physics_functions_module.pthresh(
            physics_variables.dene,
            physics_variables.dnla,
            physics_variables.bt,
            physics_variables.rmajor,
            physics_variables.kappa,
            physics_variables.sarea,
            physics_variables.aion,
            physics_variables.aspect,
        )

        #  Enforced L-H power threshold value (if constraint 15 is turned on)

        physics_variables.plhthresh = physics_variables.pthrmw[
            physics_variables.ilhthresh - 1
        ]

        #  Power transported to the divertor by charged particles,
        #  i.e. excludes neutrons and radiation, and also NBI orbit loss power,
        #  which is assumed to be absorbed by the first wall
        if physics_variables.ignite == 0:
            pinj = current_drive_variables.pinjmw
        else:
            pinj = 0.0e0

        physics_variables.pdivt = (
            physics_variables.falpha * physics_variables.palpmw
            + physics_variables.pchargemw
            + pinj
            + physics_variables.pohmmw
            - physics_variables.pradmw
        )

        #  The following line is unphysical, but prevents -ve sqrt argument
        #  Should be obsolete if constraint eqn 17 is turned on
        physics_variables.pdivt = max(0.001e0, physics_variables.pdivt)

        # if double null configuration share the power
        # over the upper and lower divertor, where physics_variables.ftar gives
        # the factor of power conducted to the lower divertor
        if physics_variables.idivrt == 2:
            physics_variables.pdivl = physics_variables.ftar * physics_variables.pdivt
            physics_variables.pdivu = (
                1.0e0 - physics_variables.ftar
            ) * physics_variables.pdivt
            physics_variables.pdivmax = max(
                physics_variables.pdivl, physics_variables.pdivu
            )

        # Resistive diffusion time = current penetration time ~ mu0.a^2/resistivity
        physics_variables.res_time = physics_functions_module.res_diff_time()

        #  Power transported to the first wall by escaped alpha particles
        physics_variables.palpfwmw = physics_variables.palpmw * (
            1.0e0 - physics_variables.falpha
        )

        #  Density limit
        physics_variables.dlimit, physics_variables.dnelimt = physics_module.culdlm(
            physics_variables.bt,
            physics_variables.idensl,
            physics_variables.pdivt,
            physics_variables.plascur,
            divertor_variables.prn1,
            physics_variables.qstar,
            physics_variables.q95,
            physics_variables.rmajor,
            physics_variables.rminor,
            physics_variables.sarea,
            physics_variables.zeff,
        )

        #  Calculate transport losses and energy confinement time using the
        #  chosen scaling law
        (
            physics_variables.kappaa,
            physics_variables.ptrepv,
            physics_variables.ptripv,
            physics_variables.tauee,
            physics_variables.taueff,
            physics_variables.tauei,
            physics_variables.powerht,
        ) = physics_module.pcond(
            physics_variables.afuel,
            physics_variables.palpmw,
            physics_variables.aspect,
            physics_variables.bt,
            physics_variables.dnitot,
            physics_variables.dene,
            physics_variables.dnla,
            physics_variables.eps,
            physics_variables.hfact,
            physics_variables.iinvqd,
            physics_variables.isc,
            physics_variables.ignite,
            physics_variables.kappa,
            physics_variables.kappa95,
            physics_variables.pchargemw,
            current_drive_variables.pinjmw,
            physics_variables.plascur,
            physics_variables.pcoreradpv,
            physics_variables.rmajor,
            physics_variables.rminor,
            physics_variables.te,
            physics_variables.ten,
            physics_variables.tin,
            physics_variables.q95,
            physics_variables.qstar,
            physics_variables.vol,
            physics_variables.xarea,
            physics_variables.zeff,
        )

        physics_variables.ptremw = physics_variables.ptrepv * physics_variables.vol
        physics_variables.ptrimw = physics_variables.ptripv * physics_variables.vol
        #  Total transport power from scaling law (MW)
        # pscalingmw = physics_variables.ptremw + physics_variables.ptrimw #KE - why is this commented?

        # Calculate Volt-second requirements
        (
            physics_variables.phiint,
            physics_variables.rlp,
            physics_variables.vsbrn,
            physics_variables.vsind,
            physics_variables.vsres,
            physics_variables.vsstt,
        ) = physics_module.vscalc(
            physics_variables.csawth,
            physics_variables.eps,
            physics_variables.facoh,
            physics_variables.gamma,
            physics_variables.kappa,
            physics_variables.rmajor,
            physics_variables.rplas,
            physics_variables.plascur,
            times_variables.theat,
            times_variables.tburn,
            physics_variables.rli,
        )

        #  Calculate auxiliary physics related information
        sbar = 1.0e0
        (
            physics_variables.burnup,
            physics_variables.dntau,
            physics_variables.figmer,
            physics_module.fusrat,
            physics_variables.qfuel,
            physics_variables.rndfuel,
            physics_variables.taup,
        ) = physics_module.phyaux(
            physics_variables.aspect,
            physics_variables.dene,
            physics_variables.deni,
            physics_variables.fusionrate,
            physics_variables.alpharate,
            physics_variables.plascur,
            sbar,
            physics_variables.dnalp,
            physics_variables.taueff,
            physics_variables.vol,
        )

        # ptremw = physics_variables.ptrepv*physics_variables.vol
        # ptrimw = physics_variables.ptripv*physics_variables.vol
        #  Total transport power from scaling law (MW)
        physics_variables.pscalingmw = (
            physics_variables.ptremw + physics_variables.ptrimw
        )

        #  Calculate physics_variables.beta limit
        if physics_variables.iprofile == 0:
            if physics_variables.gtscale == 1:
                #  Original scaling law
                physics_variables.dnbeta = 2.7e0 * (
                    1.0e0 + 5.0e0 * physics_variables.eps**3.5e0
                )

            if physics_variables.gtscale == 2:
                # See Issue #1439
                # physics_variables.dnbeta found from physics_variables.aspect ratio scaling on p32 of Menard:
                # Menard, et al. "Fusion Nuclear Science Facilities
                # and Pilot Plants Based on the Spherical Tokamak."
                # Nucl. Fusion, 2016, 44.
                physics_variables.dnbeta = (
                    3.12e0 + 3.5e0 * physics_variables.eps**1.7e0
                )

        else:
            #  Relation between physics_variables.beta limit and plasma internal inductance
            #  Hartmann and Zohm
            physics_variables.dnbeta = 4.0e0 * physics_variables.rli

        physics_variables.betalim = physics_module.culblm(
            physics_variables.bt,
            physics_variables.dnbeta,
            physics_variables.plascur,
            physics_variables.rminor,
        )

        # MDK
        #  Nominal mean photon wall load on entire first wall area including divertor and beam holes
        #  Note that 'fwarea' excludes these, so they have been added back in.
        if physics_variables.iwalld == 1:
            physics_variables.photon_wall = (
                physics_variables.ffwal
                * physics_variables.pradmw
                / physics_variables.sarea
            )
        else:
            if physics_variables.idivrt == 2:
                # Double Null configuration in - including SoL radiation
                physics_variables.photon_wall = (
                    1.0e0 - fwbs_variables.fhcd - 2.0e0 * fwbs_variables.fdiv
                ) * physics_variables.pradmw / build_variables.fwarea + (
                    1.0e0 - fwbs_variables.fhcd - 2.0e0 * fwbs_variables.fdiv
                ) * physics_variables.rad_fraction_sol * physics_variables.pdivt / (
                    build_variables.fwarea
                )
            else:
                # Single null configuration - including SoL radaition
                physics_variables.photon_wall = (
                    1.0e0 - fwbs_variables.fhcd - fwbs_variables.fdiv
                ) * physics_variables.pradmw / build_variables.fwarea + (
                    1.0e0 - fwbs_variables.fhcd - fwbs_variables.fdiv
                ) * physics_variables.rad_fraction_sol * physics_variables.pdivt / build_variables.fwarea

        constraint_variables.peakradwallload = (
            physics_variables.photon_wall * constraint_variables.peakfactrad
        )

        # Calculate the target imbalances
        # find the total power into the targets
        physics_module.ptarmw = physics_variables.pdivt * (
            1.0e0 - physics_variables.rad_fraction_sol
        )
        # use physics_variables.ftar to find deltarsep
        # Parameters taken from double null machine
        # D. Brunner et al
        physics_module.lambdaio = 1.57e-3

        # Issue #1559 Infinities in physics_module.drsep when running single null in a double null machine
        # C W Ashe
        if physics_variables.ftar < 4.5e-5:
            physics_module.drsep = 1.5e-2
        elif physics_variables.ftar > (1.0e0 - 4.5e-5):
            physics_module.drsep = -1.5e-2
        else:
            physics_module.drsep = (
                -2.0e0 * 1.5e-3 * math.atanh(2.0e0 * (physics_variables.ftar - 0.5e0))
            )
        # Model Taken from D3-D paper for conventional divertor
        # Journal of Nuclear Materials
        # Volumes 290–293, March 2001, Pages 935-939
        # Find the innner and outer lower target imbalance

        physics_module.fio = 0.16e0 + (0.16e0 - 0.41e0) * (
            1.0e0
            - (
                2.0e0
                / (
                    1.0e0
                    + numpy.exp(
                        -((physics_module.drsep / physics_module.lambdaio) ** 2)
                    )
                )
            )
        )
        if physics_variables.idivrt == 2:
            # Double Null configuration
            # Find all the power fractions accross the targets
            # Taken from D3-D conventional divertor design
            physics_module.fli = physics_variables.ftar * physics_module.fio
            physics_module.flo = physics_variables.ftar * (1.0e0 - physics_module.fio)
            physics_module.fui = (1.0e0 - physics_variables.ftar) * physics_module.fio
            physics_module.fuo = (1.0e0 - physics_variables.ftar) * (
                1.0e0 - physics_module.fio
            )
            # power into each target
            physics_module.plimw = physics_module.fli * physics_module.ptarmw
            physics_module.plomw = physics_module.flo * physics_module.ptarmw
            physics_module.puimw = physics_module.fui * physics_module.ptarmw
            physics_module.puomw = physics_module.fuo * physics_module.ptarmw
        else:
            # Single null configuration
            physics_module.fli = physics_module.fio
            physics_module.flo = 1.0e0 - physics_module.fio
            # power into each target
            physics_module.plimw = physics_module.fli * physics_module.ptarmw
            physics_module.plomw = physics_module.flo * physics_module.ptarmw

        # Calculate some derived quantities that may not have been defined earlier
        physics_module.total_loss_power = 1e6 * (
            physics_variables.falpha * physics_variables.palpmw
            + physics_variables.pchargemw
            + physics_variables.pohmmw
            + current_drive_variables.pinjmw
        )
        physics_module.rad_fraction_lcfs = (
            1.0e6 * physics_variables.pradmw / physics_module.total_loss_power
        )
        physics_variables.rad_fraction_total = (
            physics_module.rad_fraction_lcfs
            + (1.0e0 - physics_module.rad_fraction_lcfs)
            * physics_variables.rad_fraction_sol
        )
        physics_variables.pradsolmw = (
            physics_variables.rad_fraction_sol * physics_variables.pdivt
        )
        physics_module.total_plasma_internal_energy = (
            1.5e0
            * physics_variables.beta
            * physics_variables.btot
            * physics_variables.btot
            / (2.0e0 * constants.rmu0)
            * physics_variables.vol
        )

        physics_module.total_energy_conf_time = (
            physics_module.total_plasma_internal_energy
            / physics_module.total_loss_power
        )

        if any(numerics.icc == 78):
            po.write(
                self.outfile,
                (
                    f"reinke t and fz, physics = {physics_variables.tesep} , {reinke_variables.fzmin}"
                ),
            )
            # fsep = physics_variables.nesep / physics_variables.dene
            fgw = physics_variables.dlimit(7) / physics_variables.dene
            # calculate separatrix temperature, if Reinke criterion is used
            physics_variables.tesep = reinke_module.reinke_tsep(
                physics_variables.bt,
                constraint_variables.flhthresh,
                physics_variables.q95,
                physics_variables.rmajor,
                physics_variables.eps,
                fgw,
                physics_variables.kappa,
                reinke_variables.lhat,
            )

            if reinke_variables.fzmin >= 1.0e0:
                error_handling.report_error(217)

            po.write(
                self.outfile,
                (
                    f" 'fzactual, frac, reinke_variables.impvardiv = {reinke_variables.fzactual}, {impurity_radiation_module.impurity_arr_frac(reinke_variables.impvardiv)}, {reinke_variables.impvardiv}"
                ),
            )

    def bootstrap_fraction_iter89(
        self, aspect, beta, bt, cboot, plascur, q95, q0, rmajor, vol
    ):
        """Original ITER calculation of bootstrap-driven fraction
        of the plasma current.
        author: P J Knight, CCFE, Culham Science Centre
        aspect  : input real : plasma aspect ratio
        beta    : input real : plasma total beta
        bt      : input real : toroidal field on axis (T)
        cboot   : input real : bootstrap current fraction multiplier
        plascur : input real : plasma current (A)
        q95     : input real : safety factor at 95% surface
        q0      : input real : central safety factor
        rmajor  : input real : plasma major radius (m)
        vol     : input real : plasma volume (m3)
        This routine performs the original ITER calculation of the
        plasma current bootstrap fraction.
        ITER Physics Design Guidelines: 1989 [IPDG89], N. A. Uckan et al,
        ITER Documentation Series No.10, IAEA/ITER/DS/10, IAEA, Vienna, 1990
        """
        xbs = min(10, q95 / q0)
        cbs = cboot * (1.32 - 0.235 * xbs + 0.0185 * xbs**2)
        bpbs = (
            constants.rmu0
            * plascur
            / (2 * numpy.pi * numpy.sqrt(vol / (2 * numpy.pi**2 * rmajor)))
        )
        betapbs = beta * bt**2 / bpbs**2

        if betapbs <= 0.0:  # only possible if beta <= 0.0
            return 0.0
        return cbs * (betapbs / numpy.sqrt(aspect)) ** 1.3

    def bootstrap_fraction_nevins(
        self,
        alphan,
        alphat,
        betat,
        bt,
        dene,
        plascur,
        q95,
        q0,
        rmajor,
        rminor,
        ten,
        zeff,
    ):
        """Bootstrap current fraction from Nevins et al scaling
        author: P J Knight, CCFE, Culham Science Centre
        alphan : input real :  density profile index
        alphat : input real :  temperature profile index
        betat  : input real :  total plasma beta (with respect to the toroidal
        field)
        bt     : input real :  toroidal field on axis (T)
        dene   : input real :  electron density (/m3)
        plascur: input real :  plasma current (A)
        q0     : input real :  central safety factor
        q95    : input real :  safety factor at 95% surface
        rmajor : input real :  plasma major radius (m)
        rminor : input real :  plasma minor radius (m)
        ten    : input real :  density weighted average plasma temperature (keV)
        zeff   : input real :  plasma effective charge
        This function calculates the bootstrap current fraction,
        using the Nevins et al method, 4/11/90.
        AEA FUS 251: A User's Guide to the PROCESS Systems Code
        """
        # Calculate peak electron beta

        betae0 = (
            physics_variables.ne0
            * physics_variables.te0
            * 1.0e3
            * constants.echarge
            / (bt**2 / (2.0 * constants.rmu0))
        )

        # Call integration routine

        ainteg, _ = integrate.quad(
            lambda y: self.bsinteg(
                y, dene, ten, bt, rminor, rmajor, zeff, alphat, alphan, q0, q95, betat
            ),
            0,
            0.999,
            epsabs=0.001,
            epsrel=0.001,
        )

        # Calculate bootstrap current and fraction

        aibs = 2.5 * betae0 * rmajor * bt * q95 * ainteg
        return 1.0e6 * aibs / plascur

    def bsinteg(
        self, y, dene, ten, bt, rminor, rmajor, zeff, alphat, alphan, q0, q95, betat
    ):
        """Integrand function for Nevins et al bootstrap current scaling
        author: P J Knight, CCFE, Culham Science Centre
        y : input real : abscissa of integration, = normalised minor radius
        This function calculates the integrand function for the
        Nevins et al bootstrap current scaling, 4/11/90.
        AEA FUS 251: A User's Guide to the PROCESS Systems Code
        """
        #  Constants for fit to q-profile

        c1 = 1.0
        c2 = 1.0
        c3 = 1.0

        #  Compute average electron beta

        betae = (
            dene * ten * 1.0e3 * constants.echarge / (bt**2 / (2.0 * constants.rmu0))
        )

        nabla = rminor * numpy.sqrt(y) / rmajor
        x = (1.46 * numpy.sqrt(nabla) + 2.4 * nabla) / (1.0 - nabla) ** 1.5
        z = zeff
        d = (
            1.414 * z
            + z * z
            + x * (0.754 + 2.657 * z + 2.0 * z * z)
            + x * x * (0.348 + 1.243 * z + z * z)
        )
        al2 = -x * (0.884 + 2.074 * z) / d
        a2 = alphat * (1.0 - y) ** (alphan + alphat - 1.0)
        alphai = -1.172 / (1.0 + 0.462 * x)
        a1 = (alphan + alphat) * (1.0 - y) ** (alphan + alphat - 1.0)
        al1 = x * (0.754 + 2.21 * z + z * z + x * (0.348 + 1.243 * z + z * z)) / d

        #  q-profile

        q = q0 + (q95 - q0) * (c1 * y + c2 * y * y + c3 * y**3) / (c1 + c2 + c3)

        pratio = (betat - betae) / betae

        return (q / q95) * (al1 * (a1 + pratio * (a1 + alphai * a2)) + al2 * a2)

    def bootstrap_fraction_wilson(
        self, alphaj, alphap, alphat, betpth, q0, qpsi, rmajor, rminor
    ):
        """Bootstrap current fraction from Wilson et al scaling
        author: P J Knight, CCFE, Culham Science Centre
        alphaj  : input real :  current profile index
        alphap  : input real :  pressure profile index
        alphat  : input real :  temperature profile index
        beta    : input real :  total beta
        betpth  : input real :  thermal component of poloidal beta
        q0      : input real :  safety factor on axis
        qpsi    : input real :  edge safety factor
        rmajor  : input real :  major radius (m)
        rminor  : input real :  minor radius (m)
        This function calculates the bootstrap current fraction
        using the numerically fitted algorithm written by Howard Wilson.
        AEA FUS 172: Physics Assessment for the European Reactor Study
        H. R. Wilson, Nuclear Fusion <B>32</B> (1992) 257
        """
        term1 = numpy.log(0.5)
        term2 = numpy.log(q0 / qpsi)

        termp = 1.0 - 0.5 ** (1.0 / alphap)
        termt = 1.0 - 0.5 ** (1.0 / alphat)
        termj = 1.0 - 0.5 ** (1.0 / alphaj)

        alfpnw = term1 / numpy.log(numpy.log((q0 + (qpsi - q0) * termp) / qpsi) / term2)
        alftnw = term1 / numpy.log(numpy.log((q0 + (qpsi - q0) * termt) / qpsi) / term2)
        aj = term1 / numpy.log(numpy.log((q0 + (qpsi - q0) * termj) / qpsi) / term2)

        #  Crude check for NaN errors or other illegal values...

        if numpy.isnan(aj) or numpy.isnan(alfpnw) or numpy.isnan(alftnw) or aj < 0:
            error_handling.fdiags[0] = aj
            error_handling.fdiags[1] = alfpnw
            error_handling.fdiags[2] = alftnw
            error_handling.fdiags[3] = aj

            error_handling.report_error(76)

        #  Ratio of ionic charge to electron charge

        z = 1.0

        #  Inverse aspect ratio: r2 = maximum plasma radius, r1 = minimum

        r2 = rmajor + rminor
        r1 = rmajor - rminor
        eps1 = (r2 - r1) / (r2 + r1)

        #  Coefficients fitted using least squares techniques

        saj = numpy.sqrt(aj)

        a = numpy.array(
            [
                1.41 * (1.0 - 0.28 * saj) * (1.0 + 0.12 / z),
                0.36 * (1.0 - 0.59 * saj) * (1.0 + 0.8 / z),
                -0.27 * (1.0 - 0.47 * saj) * (1.0 + 3.0 / z),
                0.0053 * (1.0 + 5.0 / z),
                -0.93 * (1.0 - 0.34 * saj) * (1.0 + 0.15 / z),
                -0.26 * (1.0 - 0.57 * saj) * (1.0 - 0.27 * z),
                0.064 * (1.0 - 0.6 * aj + 0.15 * aj * aj) * (1.0 + 7.6 / z),
                -0.0011 * (1.0 + 9.0 / z),
                -0.33 * (1.0 - aj + 0.33 * aj * aj),
                -0.26 * (1.0 - 0.87 / saj - 0.16 * aj),
                -0.14 * (1.0 - 1.14 / saj - 0.45 * saj),
                -0.0069,
            ]
        )

        seps1 = numpy.sqrt(eps1)

        b = numpy.array(
            [
                1.0,
                alfpnw,
                alftnw,
                alfpnw * alftnw,
                seps1,
                alfpnw * seps1,
                alftnw * seps1,
                alfpnw * alftnw * seps1,
                eps1,
                alfpnw * eps1,
                alftnw * eps1,
                alfpnw * alftnw * eps1,
            ]
        )

        #  Empirical bootstrap current fraction

        return seps1 * betpth * (a * b).sum()

    def diamagnetic_fraction_hender(self, beta):
        """author: S.I. Muldrew, CCFE, Culham Science Centre
        Diamagnetic contribution at tight aspect ratio.
        Tim Hender fit
        """
        return beta / 2.8

    def diamagnetic_fraction_scene(self, beta, q95, q0):
        """author: S.I. Muldrew, CCFE, Culham Science Centre
        Diamagnetic fraction based on SCENE fit by Tim Hender
        See Issue #992
        """
        return beta * (0.1 * q95 / q0 + 0.44) * 4.14e-1

    def ps_fraction_scene(self, beta):
        """author: S.I. Muldrew, CCFE, Culham Science Centre
        Pfirsch-Schlüter fraction based on SCENE fit by Tim Hender
        See Issue #992
        """
        return -9e-2 * beta

    def bootstrap_fraction_sauter(self):
        """Bootstrap current fraction from Sauter et al scaling
        author: P J Knight, CCFE, Culham Science Centre
        None
        This function calculates the bootstrap current fraction
        using the Sauter, Angioni and Lin-Liu scaling.
        <P>The code was supplied by Emiliano Fable, IPP Garching
        (private communication).
        O. Sauter, C. Angioni and Y. R. Lin-Liu,
        Physics of Plasmas <B>6</B> (1999) 2834
        O. Sauter, C. Angioni and Y. R. Lin-Liu, (ERRATA)
        Physics of Plasmas <B>9</B> (2002) 5140
        """
        NR = 200

        roa = numpy.arange(1, NR + 1, step=1) / NR

        rho = numpy.sqrt(physics_variables.xarea / numpy.pi) * roa
        sqeps = numpy.sqrt(roa * (physics_variables.rminor / physics_variables.rmajor))

        ne = 1e-19 * numpy.vectorize(
            lambda r: profiles_module.nprofile(
                r,
                physics_variables.rhopedn,
                physics_variables.ne0,
                physics_variables.neped,
                physics_variables.nesep,
                physics_variables.alphan,
            )
        )(roa)
        ni = (physics_variables.dnitot / physics_variables.dene) * ne
        tempe = numpy.vectorize(
            lambda r: profiles_module.tprofile(
                r,
                physics_variables.rhopedt,
                physics_variables.te0,
                physics_variables.teped,
                physics_variables.tesep,
                physics_variables.alphat,
                physics_variables.tbeta,
            )
        )(roa)
        tempi = (physics_variables.ti / physics_variables.te) * tempe

        zef = numpy.full_like(
            tempi, physics_variables.zeff
        )  # Flat Zeff profile assumed

        # mu = 1/safety factor
        # Parabolic q profile assumed
        mu = 1 / (
            physics_variables.q0
            + (physics_variables.q - physics_variables.q0) * roa**2
        )
        amain = numpy.full_like(mu, physics_variables.afuel)
        zmain = numpy.full_like(mu, 1.0 + physics_variables.fhe3)

        if ne[NR - 1] == 0.0:
            ne[NR - 1] = 1e-4 * ne[NR - 2]
            ni[NR - 1] = 1e-4 * ni[NR - 2]

        if tempe[NR - 1] == 0.0:
            tempe[NR - 1] = 1e-4 * tempe[NR - 2]
            tempi[NR - 1] = 1e-4 * tempi[NR - 2]

        # Calculate total bootstrap current (MA) by summing along profiles
        iboot = 0.0
        for ir in range(0, NR):
            if ir + 1 == NR:
                jboot = 0.0
                da = 0.0
            else:
                drho = rho[ir + 1] - rho[ir]
                da = 2 * numpy.pi * rho[ir] * drho  # area of annulus

                dlogte_drho = (numpy.log(tempe[ir + 1]) - numpy.log(tempe[ir])) / drho
                dlogti_drho = (numpy.log(tempi[ir + 1]) - numpy.log(tempi[ir])) / drho
                dlogne_drho = (numpy.log(ne[ir + 1]) - numpy.log(ne[ir])) / drho

                # The factor of 0.5 below arises because in ASTRA the logarithms
                # are coded as (e.g.):  (Te(j+1)-Te(j))/(Te(j+1)+Te(j)), which
                # actually corresponds to grad(log(Te))/2. So the factors dcsa etc.
                # are a factor two larger than one might otherwise expect.
                jboot = 0.5 * (
                    physics_module.dcsa(
                        ir + 1,
                        NR,
                        physics_variables.rmajor,
                        physics_variables.bt,
                        physics_variables.triang,
                        ne,
                        ni,
                        tempe,
                        tempi,
                        mu,
                        rho,
                        zef,
                        sqeps,
                    )
                    * dlogne_drho
                    + physics_module.hcsa(
                        ir + 1,
                        NR,
                        physics_variables.rmajor,
                        physics_variables.bt,
                        physics_variables.triang,
                        ne,
                        ni,
                        tempe,
                        tempi,
                        mu,
                        rho,
                        zef,
                        sqeps,
                    )
                    * dlogte_drho
                    + physics_module.xcsa(
                        ir + 1,
                        NR,
                        physics_variables.rmajor,
                        physics_variables.bt,
                        physics_variables.triang,
                        mu,
                        sqeps,
                        tempi,
                        tempe,
                        amain,
                        zmain,
                        ni,
                        ne,
                        rho,
                        zef,
                    )
                    * dlogti_drho
                )
                jboot = (
                    -physics_variables.bt
                    / (0.2 * numpy.pi * physics_variables.rmajor)
                    * rho[ir]
                    * mu[ir]
                    * jboot
                )  # MA/m2

            iboot += da * jboot
        return 1.0e6 * iboot / physics_variables.plascur

    def eped_warning(self):
        eped_warning = ""

        if (physics_variables.triang < 0.399e0) or (physics_variables.triang > 0.601e0):
            eped_warning += f"{physics_variables.triang = }"

        if (physics_variables.kappa < 1.499e0) or (physics_variables.kappa > 2.001e0):
            eped_warning += f"{physics_variables.kappa = }"

        if (physics_variables.plascur < 9.99e6) or (
            physics_variables.plascur > 20.01e6
        ):
            eped_warning += f"{physics_variables.plascur = }"

        if (physics_variables.rmajor < 6.99e0) or (physics_variables.rmajor > 11.01e0):
            eped_warning += f"{physics_variables.rmajor = }"

        if (physics_variables.rminor < 1.99e0) or (physics_variables.rminor > 3.501e0):
            eped_warning += f"{physics_variables.rminor = }"

        if (physics_variables.normalised_total_beta < 1.99e0) or (
            physics_variables.normalised_total_beta > 3.01e0
        ):
            eped_warning += f"{physics_variables.normalised_total_beta = }"

        if physics_variables.tesep > 0.5:
            eped_warning += f"{physics_variables.tesep = }"

        return eped_warning

    def outtim(self):
        po.oheadr(self.outfile, "Times")

        po.ovarrf(
            self.outfile,
            "Initial charge time for CS from zero current (s)",
            "(tramp)",
            times_variables.tramp,
        )
        po.ovarrf(
            self.outfile,
            "Plasma current ramp-up time (s)",
            "(tohs)",
            times_variables.tohs,
        )
        po.ovarrf(self.outfile, "Heating time (s)", "(theat)", times_variables.theat)
        po.ovarre(
            self.outfile, "Burn time (s)", "(tburn)", times_variables.tburn, "OP "
        )
        po.ovarrf(
            self.outfile,
            "Reset time to zero current for CS (s)",
            "(tqnch)",
            times_variables.tqnch,
        )
        po.ovarrf(
            self.outfile, "Time between pulses (s)", "(tdwell)", times_variables.tdwell
        )
        po.oblnkl(self.outfile)
        po.ovarre(
            self.outfile,
            "Total plant cycle time (s)",
            "(tcycle)",
            times_variables.tcycle,
            "OP ",
        )

    def outplas(self):
        """Subroutine to output the plasma physics information
        author: P J Knight, CCFE, Culham Science Centre
        self.outfile : input integer : Fortran output unit identifier
        This routine writes the plasma physics information
        to a file, in a tidy format.
        AEA FUS 251: A User's Guide to the PROCESS Systems Code
        """

        # ###############################################
        # Dimensionless plasma parameters. See reference below.
        physics_module.nu_star = (
            1
            / constants.rmu0
            * (15.0e0 * constants.echarge**4 * physics_variables.dlamie)
            / (4.0e0 * numpy.pi**1.5e0 * constants.epsilon0**2)
            * physics_variables.vol**2
            * physics_variables.rmajor**2
            * physics_variables.bt
            * numpy.sqrt(physics_variables.eps)
            * physics_variables.dnla**3
            * physics_variables.kappa
            / (
                physics_module.total_plasma_internal_energy**2
                * physics_variables.plascur
            )
        )

        physics_module.rho_star = numpy.sqrt(
            2.0e0
            * constants.mproton
            * physics_variables.aion
            * physics_module.total_plasma_internal_energy
            / (3.0e0 * physics_variables.vol * physics_variables.dnla)
        ) / (
            constants.echarge
            * physics_variables.bt
            * physics_variables.eps
            * physics_variables.rmajor
        )

        physics_module.beta_mcdonald = (
            4.0e0
            / 3.0e0
            * constants.rmu0
            * physics_module.total_plasma_internal_energy
            / (physics_variables.vol * physics_variables.bt**2)
        )

        po.oheadr(self.outfile, "Plasma")

        if stellarator_variables.istell == 0:
            if physics_variables.idivrt == 0:
                po.ocmmnt(self.outfile, "Plasma configuration = limiter")
            elif physics_variables.idivrt == 1:
                po.ocmmnt(self.outfile, "Plasma configuration = single null divertor")
            elif physics_variables.idivrt == 2:
                po.ocmmnt(self.outfile, "Plasma configuration = double null divertor")
            else:
                error_handling.idiags[0] = physics_variables.idivrt
                po.report_error(85)
        else:
            po.ocmmnt(self.outfile, "Plasma configuration = stellarator")

        if stellarator_variables.istell == 0:
            if physics_variables.itart == 0:
                physics_module.itart_r = physics_variables.itart
                po.ovarrf(
                    self.outfile,
                    "Tokamak aspect ratio = Conventional, itart = 0",
                    "(itart)",
                    physics_module.itart_r,
                )
            elif physics_variables.itart == 1:
                physics_module.itart_r = physics_variables.itart
                po.ovarrf(
                    self.outfile,
                    "Tokamak aspect ratio = Spherical, itart = 1",
                    "(itart)",
                    physics_module.itart_r,
                )

        po.osubhd(self.outfile, "Plasma Geometry :")
        po.ovarrf(
            self.outfile, "Major radius (m)", "(rmajor)", physics_variables.rmajor
        )
        po.ovarrf(
            self.outfile,
            "Minor radius (m)",
            "(rminor)",
            physics_variables.rminor,
            "OP ",
        )
        po.ovarrf(self.outfile, "Aspect ratio", "(aspect)", physics_variables.aspect)

        if stellarator_variables.istell == 0:
            if physics_variables.ishape in [0, 6, 8]:
                po.ovarrf(
                    self.outfile,
                    "Elongation, X-point (input value used)",
                    "(kappa)",
                    physics_variables.kappa,
                    "IP ",
                )
            elif physics_variables.ishape == 1:
                po.ovarrf(
                    self.outfile,
                    "Elongation, X-point (TART scaling)",
                    "(kappa)",
                    physics_variables.kappa,
                    "OP ",
                )
            elif physics_variables.ishape in [2, 3]:
                po.ovarrf(
                    self.outfile,
                    "Elongation, X-point (Zohm scaling)",
                    "(kappa)",
                    physics_variables.kappa,
                    "OP ",
                )
                po.ovarrf(
                    self.outfile,
                    "Zohm scaling adjustment factor",
                    "(fkzohm)",
                    physics_variables.fkzohm,
                )
            elif physics_variables.ishape in [4, 5, 7]:
                po.ovarrf(
                    self.outfile,
                    "Elongation, X-point (calculated from kappa95)",
                    "(kappa)",
                    physics_variables.kappa,
                    "OP ",
                )
            elif physics_variables.ishape == 9:
                po.ovarrf(
                    self.outfile,
                    "Elongation, X-point (calculated from aspect ratio and li(3))",
                    "(kappa)",
                    physics_variables.kappa,
                    "OP ",
                )
            elif physics_variables.ishape == 10:
                po.ovarrf(
                    self.outfile,
                    "Elongation, X-point (calculated from aspect ratio and stability margin)",
                    "(kappa)",
                    physics_variables.kappa,
                    "OP ",
                )
            elif physics_variables.ishape == 11:
                po.ovarrf(
                    self.outfile,
                    "Elongation, X-point (calculated from aspect ratio via Menard 2016)",
                    "(kappa)",
                    physics_variables.kappa,
                    "OP ",
                )
            else:
                error_handling.idiags[0] = physics_variables.ishape
                po.report_error(86)

            if physics_variables.ishape in [4, 5, 7]:
                po.ovarrf(
                    self.outfile,
                    "Elongation, 95% surface (input value used)",
                    "(kappa95)",
                    physics_variables.kappa95,
                    "IP ",
                )
            else:
                po.ovarrf(
                    self.outfile,
                    "Elongation, 95% surface (calculated from physics_variables.kappa)",
                    "(kappa95)",
                    physics_variables.kappa95,
                    "OP ",
                )

            po.ovarrf(
                self.outfile,
                "Elongation, area ratio calc.",
                "(kappaa)",
                physics_variables.kappaa,
                "OP ",
            )

            if physics_variables.ishape in [0, 2, 6, 8, 9, 10, 11]:
                po.ovarrf(
                    self.outfile,
                    "Triangularity, X-point (input value used)",
                    "(triang)",
                    physics_variables.triang,
                    "IP ",
                )
            elif physics_variables.ishape == 1:
                po.ovarrf(
                    self.outfile,
                    "Triangularity, X-point (TART scaling)",
                    "(triang)",
                    physics_variables.triang,
                    "OP ",
                )
            else:
                po.ovarrf(
                    self.outfile,
                    "Triangularity, X-point (calculated from triang95)",
                    "(triang)",
                    physics_variables.triang,
                    "OP ",
                )

            if physics_variables.ishape in [3, 4, 5, 7]:
                po.ovarrf(
                    self.outfile,
                    "Triangularity, 95% surface (input value used)",
                    "(triang95)",
                    physics_variables.triang95,
                    "IP ",
                )
            else:
                po.ovarrf(
                    self.outfile,
                    "Triangularity, 95% surface (calculated from triang)",
                    "(triang95)",
                    physics_variables.triang95,
                    "OP ",
                )

            po.ovarrf(
                self.outfile,
                "Plasma poloidal perimeter (m)",
                "(pperim)",
                physics_variables.pperim,
                "OP ",
            )

            po.ovarrf(
                self.outfile,
                "Plasma cross-sectional area (m2)",
                "(xarea)",
                physics_variables.xarea,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Plasma surface area (m2)",
                "(sarea)",
                physics_variables.sarea,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Plasma volume (m3)",
                "(vol)",
                physics_variables.vol,
                "OP ",
            )

            po.osubhd(self.outfile, "Current and Field :")

            if stellarator_variables.istell == 0:
                if physics_variables.iprofile == 0:
                    po.ocmmnt(
                        self.outfile,
                        "Consistency between q0,q,alphaj,rli,dnbeta is not enforced",
                    )
                else:
                    po.ocmmnt(
                        self.outfile,
                        "Consistency between q0,q,alphaj,rli,dnbeta is enforced",
                    )

                po.oblnkl(self.outfile)
                po.ovarin(
                    self.outfile,
                    "Plasma current scaling law used",
                    "(icurr)",
                    physics_variables.icurr,
                )

                po.ovarrf(
                    self.outfile,
                    "Plasma current (MA)",
                    "(plascur/1e6)",
                    physics_variables.plascur / 1.0e6,
                    "OP ",
                )

                if physics_variables.iprofile == 1:
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

                po.ovarrf(
                    self.outfile,
                    "Plasma internal inductance, li",
                    "(rli)",
                    physics_variables.rli,
                    "OP ",
                )
                po.ovarrf(
                    self.outfile,
                    "Vertical field at plasma (T)",
                    "(bvert)",
                    physics_variables.bvert,
                    "OP ",
                )

            po.ovarrf(
                self.outfile,
                "Vacuum toroidal field at R (T)",
                "(bt)",
                physics_variables.bt,
            )
            po.ovarrf(
                self.outfile,
                "Average poloidal field (T)",
                "(bp)",
                physics_variables.bp,
                "OP ",
            )

            po.ovarrf(
                self.outfile,
                "Total field (sqrt(bp^2 + bt^2)) (T)",
                "(btot)",
                physics_variables.btot,
                "OP ",
            )

        if stellarator_variables.istell == 0:
            po.ovarrf(
                self.outfile, "Safety factor on axis", "(q0)", physics_variables.q0
            )

            if physics_variables.icurr == 2:
                po.ovarrf(
                    self.outfile, "Mean edge safety factor", "(q)", physics_variables.q
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

            if physics_variables.ishape == 1:
                po.ovarrf(
                    self.outfile,
                    "Lower limit for edge safety factor q",
                    "(qlim)",
                    physics_variables.qlim,
                    "OP ",
                )

        else:
            po.ovarrf(
                self.outfile,
                "Rotational transform",
                "(iotabar)",
                stellarator_variables.iotabar,
            )

        po.osubhd(self.outfile, "Beta Information :")

        betath = (
            physics_variables.beta - physics_variables.betaft - physics_variables.betanb
        )
        gammaft = (physics_variables.betaft + physics_variables.betanb) / betath

        po.ovarre(self.outfile, "Total plasma beta", "(beta)", physics_variables.beta)
        po.ovarre(
            self.outfile,
            "Total poloidal beta",
            "(betap)",
            physics_variables.betap,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Total toroidal beta",
            " ",
            physics_variables.beta
            * (physics_variables.btot / physics_variables.bt) ** 2,
            "OP ",
        )
        po.ovarre(
            self.outfile, "Fast alpha beta", "(betaft)", physics_variables.betaft, "OP "
        )
        po.ovarre(
            self.outfile, "Beam ion beta", "(betanb)", physics_variables.betanb, "OP "
        )
        po.ovarre(
            self.outfile,
            "(Fast alpha + beam physics_variables.beta)/(thermal physics_variables.beta)",
            "(gammaft)",
            gammaft,
            "OP ",
        )

        po.ovarre(self.outfile, "Thermal beta", " ", betath, "OP ")
        po.ovarre(
            self.outfile,
            "Thermal poloidal beta",
            " ",
            betath * (physics_variables.btot / physics_variables.bp) ** 2,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Thermal toroidal physics_variables.beta (= physics_variables.beta-exp)",
            " ",
            betath * (physics_variables.btot / physics_variables.bt) ** 2,
            "OP ",
        )

        po.ovarrf(
            self.outfile,
            "2nd stability physics_variables.beta : beta_p / (R/a)",
            "(physics_variables.eps*betap)",
            physics_variables.eps * physics_variables.betap,
            "OP ",
        )
        po.ovarrf(
            self.outfile,
            "2nd stability physics_variables.beta upper limit",
            "(epbetmax)",
            physics_variables.epbetmax,
        )

        if stellarator_variables.istell == 0:
            if physics_variables.iprofile == 1:
                po.ovarrf(
                    self.outfile,
                    "Beta g coefficient",
                    "(dnbeta)",
                    physics_variables.dnbeta,
                    "OP ",
                )
            else:
                po.ovarrf(
                    self.outfile,
                    "Beta g coefficient",
                    "(dnbeta)",
                    physics_variables.dnbeta,
                )

            po.ovarrf(
                self.outfile,
                "Normalised thermal beta",
                " ",
                1.0e8
                * betath
                * physics_variables.rminor
                * physics_variables.bt
                / physics_variables.plascur,
                "OP ",
            )

            po.ovarrf(
                self.outfile,
                "Normalised total beta",
                " ",
                physics_variables.normalised_total_beta,
                "OP ",
            )

            normalised_toroidal_beta = (
                physics_variables.normalised_total_beta
                * (physics_variables.btot / physics_variables.bt) ** 2
            )
            po.ovarrf(
                self.outfile,
                "Normalised toroidal beta",
                "(normalised_toroidal_beta)",
                normalised_toroidal_beta,
                "OP ",
            )

        if physics_variables.iculbl == 0:
            po.ovarrf(
                self.outfile,
                "Limit on total beta",
                "(betalim)",
                physics_variables.betalim,
                "OP ",
            )
        elif physics_variables.iculbl == 1:
            po.ovarrf(
                self.outfile,
                "Limit on thermal beta",
                "(betalim)",
                physics_variables.betalim,
                "OP ",
            )
        else:
            po.ovarrf(
                self.outfile,
                "Limit on thermal + NB beta",
                "(betalim)",
                physics_variables.betalim,
                "OP ",
            )

        po.ovarre(
            self.outfile,
            "Plasma thermal energy (J)",
            " ",
            1.5e0
            * betath
            * physics_variables.btot
            * physics_variables.btot
            / (2.0e0 * constants.rmu0)
            * physics_variables.vol,
            "OP ",
        )

        po.ovarre(
            self.outfile,
            "Total plasma internal energy (J)",
            "(total_plasma_internal_energy)",
            physics_module.total_plasma_internal_energy,
            "OP ",
        )

        po.osubhd(self.outfile, "Temperature and Density (volume averaged) :")
        po.ovarrf(
            self.outfile, "Electron temperature (keV)", "(te)", physics_variables.te
        )
        po.ovarrf(
            self.outfile,
            "Electron temperature on axis (keV)",
            "(te0)",
            physics_variables.te0,
            "OP ",
        )
        po.ovarrf(self.outfile, "Ion temperature (keV)", "(ti)", physics_variables.ti)
        po.ovarrf(
            self.outfile,
            "Ion temperature on axis (keV)",
            "(ti0)",
            physics_variables.ti0,
            "OP ",
        )
        po.ovarrf(
            self.outfile,
            "Electron temp., density weighted (keV)",
            "(ten)",
            physics_variables.ten,
            "OP ",
        )
        po.ovarre(
            self.outfile, "Electron density (/m3)", "(dene)", physics_variables.dene
        )
        po.ovarre(
            self.outfile,
            "Electron density on axis (/m3)",
            "(ne0)",
            physics_variables.ne0,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Line-averaged electron density (/m3)",
            "(dnla)",
            physics_variables.dnla,
            "OP ",
        )

        if stellarator_variables.istell == 0:
            po.ovarre(
                self.outfile,
                "Line-averaged electron density / Greenwald density",
                "(dnla_gw)",
                physics_variables.dnla / physics_variables.dlimit[6],
                "OP ",
            )

        po.ovarre(
            self.outfile,
            "Ion density (/m3)",
            "(dnitot)",
            physics_variables.dnitot,
            "OP ",
        )
        po.ovarre(
            self.outfile, "Fuel density (/m3)", "(deni)", physics_variables.deni, "OP "
        )
        po.ovarre(
            self.outfile,
            "Total impurity density with Z > 2 (no He) (/m3)",
            "(dnz)",
            physics_variables.dnz,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Helium ion density (thermalised ions only) (/m3)",
            "(dnalp)",
            physics_variables.dnalp,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Proton density (/m3)",
            "(dnprot)",
            physics_variables.dnprot,
            "OP ",
        )
        if physics_variables.protium > 1.0e-10:
            po.ovarre(
                self.outfile,
                "Seeded protium density / electron density",
                "(protium)",
                physics_variables.protium,
            )

        po.ovarre(
            self.outfile,
            "Hot beam density (/m3)",
            "(dnbeam)",
            physics_variables.dnbeam,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Density limit from scaling (/m3)",
            "(dnelimt)",
            physics_variables.dnelimt,
            "OP ",
        )
        if (numerics.ioptimz > 0) and (numerics.active_constraints[4]):
            po.ovarre(
                self.outfile,
                "Density limit (enforced) (/m3)",
                "(numerics.boundu(9)*dnelimt)",
                numerics.boundu[8] * physics_variables.dnelimt,
                "OP ",
            )

        po.ovarre(
            self.outfile,
            "Helium ion density (thermalised ions only) / electron density",
            "(ralpne)",
            physics_variables.ralpne,
        )
        po.oblnkl(self.outfile)

        po.ocmmnt(self.outfile, "Impurities")
        po.oblnkl(self.outfile)
        po.ocmmnt(self.outfile, "Plasma ion densities / electron density:")

        for imp in range(impurity_radiation_module.nimp):
            # MDK Update fimp, as this will make the ITV output work correctly.
            impurity_radiation_module.fimp[
                imp
            ] = impurity_radiation_module.impurity_arr_frac[imp]
            str1 = (
                f2py_compatible_to_string(
                    impurity_radiation_module.impurity_arr_label[imp]
                )
                + " concentration"
            )
            str2 = f"(fimp({imp+1}))"
            # MDK Add output flag for H which is calculated.
            if imp == 0:
                po.ovarre(
                    self.outfile, str1, str2, impurity_radiation_module.fimp[imp], "OP "
                )
            else:
                po.ovarre(self.outfile, str1, str2, impurity_radiation_module.fimp[imp])

        po.ovarre(
            self.outfile,
            "Average mass of all ions (amu)",
            "(aion)",
            physics_variables.aion,
            "OP ",
        )
        # MDK Say which impurity is varied, if iteration variable fimpvar (102) is turned on
        # if (any(ixc == 102)) :
        #    call ovarst(self.outfile,'Impurity used as an iteration variable' , '', '"' // impurity_arr(impvar)%label // '"')
        #    po.ovarre(self.outfile,'Fractional density of variable impurity (ion / electron density)','(fimpvar)',fimpvar)
        #
        po.oblnkl(self.outfile)
        po.ovarrf(
            self.outfile, "Effective charge", "(zeff)", physics_variables.zeff, "OP "
        )

        # Issue #487.  No idea what zeffai is.
        # I haven't removed it as it is used in subroutine rether,
        #   (routine to find the equilibration power between the ions and electrons)
        # po.ovarrf(self.outfile,'Mass weighted effective charge','(zeffai)',zeffai, 'OP ')

        po.ovarrf(
            self.outfile, "Density profile factor", "(alphan)", physics_variables.alphan
        )
        po.ovarin(
            self.outfile,
            "Plasma profile model",
            "(ipedestal)",
            physics_variables.ipedestal,
        )

        if physics_variables.ipedestal >= 1:
            if physics_variables.ne0 < physics_variables.neped:
                error_handling.report_error(213)

            po.ocmmnt(self.outfile, "Pedestal profiles are used.")
            po.ovarrf(
                self.outfile,
                "Density pedestal r/a location",
                "(rhopedn)",
                physics_variables.rhopedn,
            )
            if physics_variables.fgwped >= 0e0:
                po.ovarre(
                    self.outfile,
                    "Electron density pedestal height (/m3)",
                    "(neped)",
                    physics_variables.neped,
                    "OP ",
                )
            else:
                po.ovarre(
                    self.outfile,
                    "Electron density pedestal height (/m3)",
                    "(neped)",
                    physics_variables.neped,
                )

            # This code is ODD# Don't change it# No explanation why fgwped and physics_variables.fgwsep
            # must be assigned to their exisiting values#
            fgwped_out = physics_variables.neped / physics_variables.dlimit[6]
            fgwsep_out = physics_variables.nesep / physics_variables.dlimit[6]
            if physics_variables.fgwped >= 0e0:
                physics_variables.fgwped = (
                    physics_variables.neped / physics_variables.dlimit[6]
                )
            if physics_variables.fgwsep >= 0e0:
                physics_variables.fgwsep = (
                    physics_variables.nesep / physics_variables.dlimit[6]
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
                "(rhopedt)",
                physics_variables.rhopedt,
            )
            # Issue #413 Pedestal scaling
            po.ovarin(
                self.outfile,
                "Pedestal scaling switch",
                "(ieped)",
                physics_variables.ieped,
            )
            if physics_variables.ieped == 1:
                po.ocmmnt(
                    self.outfile,
                    "Saarelma 6-parameter pedestal temperature scaling is ON",
                )

                if self.eped_warning() != "":
                    po.ocmmnt(
                        self.outfile,
                        "WARNING: Pedestal parameters are outside the range of applicability of the scaling:",
                    )
                    po.ocmmnt(
                        self.outfile,
                        "triang: 0.4 - 0.6; physics_variables.kappa: 1.5 - 2.0;   plascur: 10 - 20 MA, physics_variables.rmajor: 7 - 11 m;",
                    )
                    po.ocmmnt(
                        self.outfile,
                        "rminor: 2 - 3.5 m; tesep: 0 - 0.5 keV; normalised_total_beta: 2 - 3; ",
                    )
                    print(
                        "WARNING: Pedestal parameters are outside the range of applicability of the scaling:"
                    )
                    print(
                        "triang: 0.4 - 0.6; physics_variables.kappa: 1.5 - 2.0;   plascur: 10 - 20 MA, physics_variables.rmajor: 7 - 11 m;"
                    )
                    print(
                        "rminor: 2 - 3.5 m; tesep: 0 - 0.5 keV; normalised_total_beta: 2 - 3"
                    )
                    print(self.eped_warning())

            po.ovarrf(
                self.outfile,
                "Electron temp. pedestal height (keV)",
                "(teped)",
                physics_variables.teped,
            )
            if any(numerics.icc == 78):
                po.ovarrf(
                    self.outfile,
                    "Electron temp. at separatrix (keV)",
                    "(tesep)",
                    physics_variables.tesep,
                    "OP ",
                )
            else:
                po.ovarrf(
                    self.outfile,
                    "Electron temp. at separatrix (keV)",
                    "(tesep)",
                    physics_variables.tesep,
                )

            po.ovarre(
                self.outfile,
                "Electron density at separatrix (/m3)",
                "(nesep)",
                physics_variables.nesep,
            )
            po.ovarre(
                self.outfile,
                "Electron density at separatrix / nGW",
                "(fgwsep_out)",
                fgwsep_out,
            )

        # Issue 558 - addition of constraint 76 to limit the value of nesep, in proportion with the ballooning parameter and Greenwald density
        if any(numerics.icc == 76):
            po.ovarre(
                self.outfile,
                "Critical ballooning parameter value",
                "(alpha_crit)",
                physics_variables.alpha_crit,
            )
            po.ovarre(
                self.outfile,
                "Critical electron density at separatrix (/m3)",
                "(nesep_crit)",
                physics_variables.nesep_crit,
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

        if stellarator_variables.istell == 0:
            po.osubhd(self.outfile, "Density Limit using different models :")
            po.ovarre(
                self.outfile,
                "Old ASDEX model",
                "(dlimit(1))",
                physics_variables.dlimit[0],
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Borrass ITER model I",
                "(dlimit(2))",
                physics_variables.dlimit[1],
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Borrass ITER model II",
                "(dlimit(3))",
                physics_variables.dlimit[2],
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "JET edge radiation model",
                "(dlimit(4))",
                physics_variables.dlimit[3],
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "JET simplified model",
                "(dlimit(5))",
                physics_variables.dlimit[4],
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Hugill-Murakami Mq model",
                "(dlimit(6))",
                physics_variables.dlimit[5],
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Greenwald model",
                "(dlimit(7))",
                physics_variables.dlimit[6],
                "OP ",
            )

        po.osubhd(self.outfile, "Fuel Constituents :")
        po.ovarrf(
            self.outfile, "Deuterium fuel fraction", "(fdeut)", physics_variables.fdeut
        )
        po.ovarrf(
            self.outfile, "Tritium fuel fraction", "(ftrit)", physics_variables.ftrit
        )
        if physics_variables.fhe3 > 1.0e-3:
            po.ovarrf(
                self.outfile, "3-Helium fuel fraction", "(fhe3)", physics_variables.fhe3
            )

        po.osubhd(self.outfile, "Fusion Power :")
        po.ovarre(
            self.outfile,
            "Total fusion power (MW)",
            "(powfmw)",
            physics_variables.powfmw,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            " =    D-T fusion power (MW)",
            "(pdt)",
            physics_variables.pdt,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "  +   D-D fusion power (MW)",
            "(pdd)",
            physics_variables.pdd,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "  + D-He3 fusion power (MW)",
            "(pdhe3)",
            physics_variables.pdhe3,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Alpha power: total (MW)",
            "(palpmw)",
            physics_variables.palpmw,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Alpha power: beam-plasma (MW)",
            "(palpnb)",
            physics_variables.palpnb,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Neutron power (MW)",
            "(pneutmw)",
            physics_variables.pneutmw,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Charged particle power (excluding alphas) (MW)",
            "(pchargemw)",
            physics_variables.pchargemw,
            "OP ",
        )
        tot_power_plasma = (
            physics_variables.falpha * physics_variables.palpmw
            + physics_variables.pchargemw
            + physics_variables.pohmmw
            + current_drive_variables.pinjmw
        )
        po.ovarre(
            self.outfile,
            "Total power deposited in plasma (MW)",
            "(tot_power_plasma)",
            tot_power_plasma,
            "OP ",
        )
        # po.ovarre(self.outfile,'Total power deposited in plasma (MW)','()',falpha*palpmw+pchargemw+pohmmw+pinjmw, 'OP ')

        po.osubhd(self.outfile, "Radiation Power (excluding SOL):")
        po.ovarre(
            self.outfile,
            "Bremsstrahlung radiation power (MW)",
            "(pbrempv*vol)",
            physics_variables.pbrempv * physics_variables.vol,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Line radiation power (MW)",
            "(plinepv*vol)",
            physics_variables.plinepv * physics_variables.vol,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Synchrotron radiation power (MW)",
            "(psyncpv*vol)",
            physics_variables.psyncpv * physics_variables.vol,
            "OP ",
        )
        po.ovarrf(
            self.outfile,
            "Synchrotron wall reflectivity factor",
            "(ssync)",
            physics_variables.ssync,
        )
        po.ovarre(
            self.outfile,
            "Normalised minor radius defining 'core'",
            "(coreradius)",
            impurity_radiation_module.coreradius,
        )
        po.ovarre(
            self.outfile,
            "Fraction of core radiation subtracted from P_L",
            "(coreradiationfraction)",
            impurity_radiation_module.coreradiationfraction,
        )
        po.ovarre(
            self.outfile,
            "Radiation power from inner zone (MW)",
            "(pinnerzoneradmw)",
            physics_variables.pinnerzoneradmw,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Radiation power from outer zone (MW)",
            "(pouterzoneradmw)",
            physics_variables.pouterzoneradmw,
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
            "Total radiation power from inside LCFS (MW)",
            "(pradmw)",
            physics_variables.pradmw,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "LCFS radiation fraction = total radiation in LCFS / total power deposited in plasma",
            "(rad_fraction_LCFS)",
            physics_module.rad_fraction_lcfs,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Nominal mean radiation load on inside surface of reactor (MW/m2)",
            "(photon_wall)",
            physics_variables.photon_wall,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Peaking factor for radiation wall load",
            "(peakfactrad)",
            constraint_variables.peakfactrad,
            "IP ",
        )
        po.ovarre(
            self.outfile,
            "Maximum permitted radiation wall load (MW/m^2)",
            "(maxradwallload)",
            constraint_variables.maxradwallload,
            "IP ",
        )
        po.ovarre(
            self.outfile,
            "Peak radiation wall load (MW/m^2)",
            "(peakradwallload)",
            constraint_variables.peakradwallload,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Fast alpha particle power incident on the first wall (MW)",
            "(palpfwmw)",
            physics_variables.palpfwmw,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Nominal mean neutron load on inside surface of reactor (MW/m2)",
            "(wallmw)",
            physics_variables.wallmw,
            "OP ",
        )

        if stellarator_variables.istell == 0:
            po.oblnkl(self.outfile)
            po.ovarre(
                self.outfile,
                "Power incident on the divertor targets (MW)",
                "(ptarmw)",
                physics_module.ptarmw,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Fraction of power to the lower divertor",
                "(ftar)",
                physics_variables.ftar,
                "IP ",
            )
            po.ovarre(
                self.outfile,
                "Outboard side heat flux decay length (m)",
                "(lambdaio)",
                physics_module.lambdaio,
                "OP ",
            )
            if physics_variables.idivrt == 2:
                po.ovarre(
                    self.outfile,
                    "Midplane seperation of the two magnetic closed flux surfaces (m)",
                    "(drsep)",
                    physics_module.drsep,
                    "OP ",
                )

            po.ovarre(
                self.outfile,
                "Fraction of power on the inner targets",
                "(fio)",
                physics_module.fio,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Fraction of power incident on the lower inner target",
                "(fLI)",
                physics_module.fli,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Fraction of power incident on the lower outer target",
                "(fLO)",
                physics_module.flo,
                "OP ",
            )
            if physics_variables.idivrt == 2:
                po.ovarre(
                    self.outfile,
                    "Fraction of power incident on the upper inner target",
                    "(fUI)",
                    physics_module.fui,
                    "OP ",
                )
                po.ovarre(
                    self.outfile,
                    "Fraction of power incident on the upper outer target",
                    "(fUO)",
                    physics_module.fuo,
                    "OP ",
                )

            po.ovarre(
                self.outfile,
                "Power incident on the lower inner target (MW)",
                "(pLImw)",
                physics_module.plimw,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Power incident on the lower outer target (MW)",
                "(pLOmw)",
                physics_module.plomw,
                "OP ",
            )
            if physics_variables.idivrt == 2:
                po.ovarre(
                    self.outfile,
                    "Power incident on the upper innner target (MW)",
                    "(pUImw)",
                    physics_module.puimw,
                    "OP ",
                )
                po.ovarre(
                    self.outfile,
                    "Power incident on the upper outer target (MW)",
                    "(pUOmw)",
                    physics_module.puomw,
                    "OP ",
                )

        po.oblnkl(self.outfile)
        po.ovarre(
            self.outfile,
            "Ohmic heating power (MW)",
            "(pohmmw)",
            physics_variables.pohmmw,
            "OP ",
        )
        po.ovarrf(
            self.outfile,
            "Fraction of alpha power deposited in plasma",
            "(falpha)",
            physics_variables.falpha,
            "OP ",
        )
        po.ovarrf(
            self.outfile,
            "Fraction of alpha power to electrons",
            "(falpe)",
            physics_variables.falpe,
            "OP ",
        )
        po.ovarrf(
            self.outfile,
            "Fraction of alpha power to ions",
            "(falpi)",
            physics_variables.falpi,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Ion transport (MW)",
            "(ptrimw)",
            physics_variables.ptrimw,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Electron transport (MW)",
            "(ptremw)",
            physics_variables.ptremw,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Injection power to ions (MW)",
            "(pinjimw)",
            current_drive_variables.pinjimw,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Injection power to electrons (MW)",
            "(pinjemw)",
            current_drive_variables.pinjemw,
            "OP ",
        )
        if physics_variables.ignite == 1:
            po.ocmmnt(self.outfile, "  (Injected power only used for start-up phase)")

        po.ovarin(
            self.outfile,
            "Ignited plasma switch (0=not ignited, 1=ignited)",
            "(ignite)",
            physics_variables.ignite,
        )

        po.oblnkl(self.outfile)
        po.ovarre(
            self.outfile,
            "Power into divertor zone via charged particles (MW)",
            "(pdivt)",
            physics_variables.pdivt,
            "OP ",
        )

        if physics_variables.pdivt <= 0.001e0:
            error_handling.fdiags[0] = physics_variables.pdivt
            error_handling.report_error(87)
            po.oblnkl(self.outfile)
            po.ocmmnt(
                self.outfile, "  BEWARE: possible problem with high radiation power"
            )
            po.ocmmnt(
                self.outfile, "          Power into divertor zone is unrealistic;"
            )
            po.ocmmnt(self.outfile, "          divertor calculations will be nonsense#")
            po.ocmmnt(
                self.outfile, "  Set constraint 17 (Radiation fraction upper limit)."
            )
            po.oblnkl(self.outfile)

        if physics_variables.idivrt == 2:
            # Double null divertor configuration
            po.ovarre(
                self.outfile,
                "Pdivt / R ratio (MW/m) (On peak divertor)",
                "(pdivmax/physics_variables.rmajor)",
                physics_variables.pdivmax / physics_variables.rmajor,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Pdivt Bt / qAR ratio (MWT/m) (On peak divertor)",
                "(pdivmaxbt/qar)",
                (
                    (physics_variables.pdivmax * physics_variables.bt)
                    / (
                        physics_variables.q95
                        * physics_variables.aspect
                        * physics_variables.rmajor
                    )
                ),
                "OP ",
            )
        else:
            # Single null divertor configuration
            po.ovarre(
                self.outfile,
                "Psep / R ratio (MW/m)",
                "(pdivt/physics_variables.rmajor)",
                physics_variables.pdivt / physics_variables.rmajor,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Psep Bt / qAR ratio (MWT/m)",
                "(pdivtbt/qar)",
                (
                    (physics_variables.pdivt * physics_variables.bt)
                    / (
                        physics_variables.q95
                        * physics_variables.aspect
                        * physics_variables.rmajor
                    )
                ),
                "OP ",
            )

        if stellarator_variables.istell == 0:
            po.osubhd(self.outfile, "H-mode Power Threshold Scalings :")

            po.ovarre(
                self.outfile,
                "ITER 1996 scaling: nominal (MW)",
                "(pthrmw(1))",
                physics_variables.pthrmw[0],
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "ITER 1996 scaling: upper bound (MW)",
                "(pthrmw(2))",
                physics_variables.pthrmw[1],
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "ITER 1996 scaling: lower bound (MW)",
                "(pthrmw(3))",
                physics_variables.pthrmw[2],
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "ITER 1997 scaling (1) (MW)",
                "(pthrmw(4))",
                physics_variables.pthrmw[3],
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "ITER 1997 scaling (2) (MW)",
                "(pthrmw(5))",
                physics_variables.pthrmw[4],
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Martin 2008 scaling: nominal (MW)",
                "(pthrmw(6))",
                physics_variables.pthrmw[5],
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Martin 2008 scaling: 95% upper bound (MW)",
                "(pthrmw(7))",
                physics_variables.pthrmw[6],
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Martin 2008 scaling: 95% lower bound (MW)",
                "(pthrmw(8))",
                physics_variables.pthrmw[7],
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Snipes 2000 scaling: nominal (MW)",
                "(pthrmw(9))",
                physics_variables.pthrmw[8],
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Snipes 2000 scaling: upper bound (MW)",
                "(pthrmw(10))",
                physics_variables.pthrmw[9],
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Snipes 2000 scaling: lower bound (MW)",
                "(pthrmw(11))",
                physics_variables.pthrmw[10],
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Snipes 2000 scaling (closed divertor): nominal (MW)",
                "(pthrmw(12))",
                physics_variables.pthrmw[11],
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Snipes 2000 scaling (closed divertor): upper bound (MW)",
                "(pthrmw(13))",
                physics_variables.pthrmw[12],
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Snipes 2000 scaling (closed divertor): lower bound (MW)",
                "(pthrmw(14))",
                physics_variables.pthrmw[13],
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Hubbard 2012 L-I threshold - nominal (MW)",
                "(pthrmw(15))",
                physics_variables.pthrmw[14],
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Hubbard 2012 L-I threshold - lower bound (MW)",
                "(pthrmw(16))",
                physics_variables.pthrmw[15],
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Hubbard 2012 L-I threshold - upper bound (MW)",
                "(pthrmw(17))",
                physics_variables.pthrmw[16],
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Hubbard 2017 L-I threshold",
                "(pthrmw(18))",
                physics_variables.pthrmw[17],
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Martin 2008 aspect ratio corrected scaling: nominal (MW)",
                "(pthrmw(19))",
                physics_variables.pthrmw[18],
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Martin 2008 aspect ratio corrected scaling: 95% upper bound (MW)",
                "(pthrmw(20))",
                physics_variables.pthrmw[19],
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Martin 2008 aspect ratio corrected scaling: 95% lower bound (MW)",
                "(pthrmw(21))",
                physics_variables.pthrmw[20],
                "OP ",
            )
            po.oblnkl(self.outfile)
            if physics_variables.ilhthresh in [9, 10, 11]:
                if (physics_variables.bt < 0.78e0) or (physics_variables.bt > 7.94e0):
                    po.ocmmnt(
                        self.outfile,
                        "(physics_variables.bt outside Snipes 2000 fitted range)",
                    )
                    error_handling.report_error(201)

                if (physics_variables.rminor < 0.15e0) or (
                    physics_variables.rminor > 1.15e0
                ):
                    po.ocmmnt(self.outfile, "(rminor outside Snipes 2000 fitted range)")
                    error_handling.report_error(202)

                if (physics_variables.rmajor < 0.55e0) or (
                    physics_variables.rmajor > 3.37e0
                ):
                    po.ocmmnt(
                        self.outfile,
                        "(physics_variables.rmajor outside Snipes 2000 fitted range)",
                    )
                    error_handling.report_error(203)

                if (physics_variables.dnla < 0.09e20) or (
                    physics_variables.dnla > 3.16e20
                ):
                    po.ocmmnt(
                        self.outfile,
                        "(physics_variables.dnla outside Snipes 2000 fitted range)",
                    )
                    error_handling.report_error(204)

                if (physics_variables.kappa < 1.0e0) or (
                    physics_variables.kappa > 2.04e0
                ):
                    po.ocmmnt(
                        self.outfile,
                        "(physics_variables.kappa outside Snipes 2000 fitted range)",
                    )
                    error_handling.report_error(205)

                if (physics_variables.triang < 0.07e0) or (
                    physics_variables.triang > 0.74e0
                ):
                    po.ocmmnt(self.outfile, "(triang outside Snipes 2000 fitted range)")
                    error_handling.report_error(206)

            po.oblnkl(self.outfile)

            if physics_variables.ilhthresh in [12, 13, 14]:
                po.ocmmnt(
                    self.outfile,
                    "(L-H threshold for closed divertor only. Limited data used in Snipes fit)",
                )
                po.oblnkl(self.outfile)
                error_handling.report_error(207)

            if (numerics.ioptimz > 0) and (numerics.active_constraints[14]):
                po.ovarre(
                    self.outfile,
                    "L-H threshold power (enforced) (MW)",
                    "(boundl(103)*plhthresh)",
                    numerics.boundl[102] * physics_variables.plhthresh,
                    "OP ",
                )
                po.ovarre(
                    self.outfile,
                    "L-H threshold power (MW)",
                    "(plhthresh)",
                    physics_variables.plhthresh,
                    "OP ",
                )
            else:
                po.ovarre(
                    self.outfile,
                    "L-H threshold power (NOT enforced) (MW)",
                    "(plhthresh)",
                    physics_variables.plhthresh,
                    "OP ",
                )

        po.osubhd(self.outfile, "Confinement :")

        if physics_variables.ignite == 1:
            po.ocmmnt(
                self.outfile,
                "Device is assumed to be ignited for the calculation of confinement time",
            )
            po.oblnkl(self.outfile)

        po.ocmmnt(
            self.outfile,
            f"Confinement scaling law: { physics_variables.tauscl[physics_variables.isc-1]}",
        )

        po.ovarrf(
            self.outfile, "Confinement H factor", "(hfact)", physics_variables.hfact
        )
        po.ovarrf(
            self.outfile,
            "Global thermal energy confinement time (s)",
            "(taueff)",
            physics_variables.taueff,
            "OP ",
        )
        po.ovarrf(
            self.outfile,
            "Ion energy confinement time (s)",
            "(tauei)",
            physics_variables.tauei,
            "OP ",
        )
        po.ovarrf(
            self.outfile,
            "Electron energy confinement time (s)",
            "(tauee)",
            physics_variables.tauee,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "n.tau = Volume-average electron density x Energy confinement time (s/m3)",
            "(dntau)",
            physics_variables.dntau,
            "OP ",
        )
        po.ocmmnt(
            self.outfile,
            "Triple product = Vol-average electron density x Vol-average         & electron temperature x Energy confinement time:",
        )
        po.ovarre(
            self.outfile,
            "Triple product  (keV s/m3)",
            "(dntau*te)",
            physics_variables.dntau * physics_variables.te,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Transport loss power assumed in scaling law (MW)",
            "(powerht)",
            physics_variables.powerht,
            "OP ",
        )
        po.ovarin(
            self.outfile,
            "Switch for radiation loss term usage in power balance",
            "(iradloss)",
            physics_variables.iradloss,
        )
        if physics_variables.iradloss == 0:
            po.ovarre(
                self.outfile,
                "Radiation power subtracted from plasma power balance (MW)",
                "",
                physics_variables.pradmw,
                "OP ",
            )
            po.ocmmnt(self.outfile, "  (Radiation correction is total radiation power)")
        elif physics_variables.iradloss == 1:
            po.ovarre(
                self.outfile,
                "Radiation power subtracted from plasma power balance (MW)",
                "",
                physics_variables.pinnerzoneradmw,
                "OP ",
            )
            po.ocmmnt(self.outfile, "  (Radiation correction is core radiation power)")
        else:
            po.ovarre(
                self.outfile,
                "Radiation power subtracted from plasma power balance (MW)",
                "",
                0.0e0,
            )
            po.ocmmnt(self.outfile, "  (No radiation correction applied)")

        po.ovarrf(
            self.outfile,
            "Alpha particle confinement time (s)",
            "(taup)",
            physics_variables.taup,
            "OP ",
        )
        # Note alpha confinement time is no longer equal to fuel particle confinement time.
        po.ovarrf(
            self.outfile,
            "Alpha particle/energy confinement time ratio",
            "(taup/taueff)",
            physics_variables.taup / physics_variables.taueff,
            "OP ",
        )
        po.ovarrf(
            self.outfile,
            "Lower limit on taup/taueff",
            "(taulimit)",
            constraint_variables.taulimit,
        )
        po.ovarrf(
            self.outfile,
            "Total energy confinement time including radiation loss (s)",
            "(total_energy_conf_time)",
            physics_module.total_energy_conf_time,
            "OP ",
        )
        po.ocmmnt(
            self.outfile,
            "  (= stored energy including fast particles / loss power including radiation",
        )

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
                physics_module.beta_mcdonald,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Normalized ion Larmor radius",
                "(rho_star)",
                physics_module.rho_star,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Normalized collisionality",
                "(nu_star)",
                physics_module.nu_star,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Volume measure of elongation",
                "(kappaa_IPB)",
                physics_variables.kappaa_ipb,
                "OP ",
            )

            po.osubhd(self.outfile, "Plasma Volt-second Requirements :")
            po.ovarre(
                self.outfile,
                "Total volt-second requirement (Wb)",
                "(vsstt)",
                physics_variables.vsstt,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Inductive volt-seconds (Wb)",
                "(vsind)",
                physics_variables.vsind,
                "OP ",
            )
            po.ovarrf(
                self.outfile, "Ejima coefficient", "(gamma)", physics_variables.gamma
            )
            po.ovarre(
                self.outfile,
                "Start-up resistive (Wb)",
                "(vsres)",
                physics_variables.vsres,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Flat-top resistive (Wb)",
                "(vsbrn)",
                physics_variables.vsbrn,
                "OP ",
            )

            po.ovarrf(
                self.outfile,
                "bootstrap current fraction multiplier",
                "(cboot)",
                current_drive_variables.cboot,
            )
            po.ovarrf(
                self.outfile,
                "Bootstrap fraction (ITER 1989)",
                "(bscf_iter89)",
                current_drive_variables.bscf_iter89,
                "OP ",
            )

            po.ovarrf(
                self.outfile,
                "Bootstrap fraction (Sauter et al)",
                "(bscf_sauter)",
                current_drive_variables.bscf_sauter,
                "OP ",
            )

            po.ovarrf(
                self.outfile,
                "Bootstrap fraction (Nevins et al)",
                "(bscf_nevins)",
                current_drive_variables.bscf_nevins,
                "OP ",
            )
            po.ovarrf(
                self.outfile,
                "Bootstrap fraction (Wilson)",
                "(bscf_wilson)",
                current_drive_variables.bscf_wilson,
                "OP ",
            )
            po.ovarrf(
                self.outfile,
                "Diamagnetic fraction (Hender)",
                "(diacf_hender)",
                current_drive_variables.diacf_hender,
                "OP ",
            )
            po.ovarrf(
                self.outfile,
                "Diamagnetic fraction (SCENE)",
                "(diacf_scene)",
                current_drive_variables.diacf_scene,
                "OP ",
            )
            po.ovarrf(
                self.outfile,
                "Pfirsch-Schlueter fraction (SCENE)",
                "(pscf_scene)",
                current_drive_variables.pscf_scene,
                "OP ",
            )
            # Error to catch if bootstap fraction limit has been enforced
            if physics_module.err242 == 1:
                error_handling.report_error(242)

            # Error to catch if self-driven current fraction limit has been enforced
            if physics_module.err243 == 1:
                error_handling.report_error(243)

            if current_drive_variables.bscfmax < 0.0e0:
                po.ocmmnt(
                    self.outfile, "  (User-specified bootstrap current fraction used)"
                )
            elif physics_variables.ibss == 1:
                po.ocmmnt(
                    self.outfile, "  (ITER 1989 bootstrap current fraction model used)"
                )
            elif physics_variables.ibss == 2:
                po.ocmmnt(
                    self.outfile,
                    "  (Nevins et al bootstrap current fraction model used)",
                )
            elif physics_variables.ibss == 3:
                po.ocmmnt(
                    self.outfile, "  (Wilson bootstrap current fraction model used)"
                )
            elif physics_variables.ibss == 4:
                po.ocmmnt(
                    self.outfile,
                    "  (Sauter et al bootstrap current fraction model used)",
                )

            if physics_variables.idia == 0:
                po.ocmmnt(
                    self.outfile, "  (Diamagnetic current fraction not calculated)"
                )
                # Error to show if diamagnetic current is above 1% but not used
                if current_drive_variables.diacf_scene > 0.01e0:
                    error_handling.report_error(244)

            elif physics_variables.idia == 1:
                po.ocmmnt(
                    self.outfile, "  (Hender diamagnetic current fraction scaling used)"
                )
            elif physics_variables.idia == 2:
                po.ocmmnt(
                    self.outfile, "  (SCENE diamagnetic current fraction scaling used)"
                )

            if physics_variables.ips == 0:
                po.ocmmnt(
                    self.outfile, "  Pfirsch-Schluter current fraction not calculated"
                )
            elif physics_variables.ips == 1:
                po.ocmmnt(
                    self.outfile,
                    "  (SCENE Pfirsch-Schluter current fraction scaling used)",
                )

            po.ovarrf(
                self.outfile,
                "Bootstrap fraction (enforced)",
                "(bootipf.)",
                current_drive_variables.bootipf,
                "OP ",
            )
            po.ovarrf(
                self.outfile,
                "Diamagnetic fraction (enforced)",
                "(diaipf.)",
                current_drive_variables.diaipf,
                "OP ",
            )
            po.ovarrf(
                self.outfile,
                "Pfirsch-Schlueter fraction (enforced)",
                "(psipf.)",
                current_drive_variables.psipf,
                "OP ",
            )

            po.ovarre(
                self.outfile,
                "Loop voltage during burn (V)",
                "(vburn)",
                physics_variables.plascur
                * physics_variables.rplas
                * physics_variables.facoh,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Plasma resistance (ohm)",
                "(rplas)",
                physics_variables.rplas,
                "OP ",
            )

            po.ovarre(
                self.outfile,
                "Resistive diffusion time (s)",
                "(res_time)",
                physics_variables.res_time,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Plasma inductance (H)",
                "(rlp)",
                physics_variables.rlp,
                "OP ",
            )
            po.ovarrf(
                self.outfile,
                "Coefficient for sawtooth effects on burn V-s requirement",
                "(csawth)",
                physics_variables.csawth,
            )

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
            "(qfuel)",
            physics_variables.qfuel,
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

        if any(numerics.icc == 78):
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
