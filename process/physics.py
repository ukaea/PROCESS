import numpy
import math
import process.physics_functions as physics_funcs
from process.fortran import constraint_variables
from process.fortran import reinke_variables
from process.fortran import reinke_module
from process.fortran import impurity_radiation_module
from process.fortran import constants
from process.fortran import physics_functions_module
from process.fortran import physics_variables
from process.fortran import physics_module
from process.fortran import pulse_variables
from process.fortran import times_variables
from process.fortran import current_drive_variables
from process.fortran import error_handling
from process.fortran import fwbs_variables
from process.fortran import build_variables
from process.fortran import divertor_variables
from process.fortran import numerics
from process.fortran import process_output as po


class Physics:
    def __init__(self, plasma_profile, current_drive):
        self.outfile = constants.nout
        self.plasma_profile = plasma_profile
        self.current_drive = current_drive

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
        current_drive_variables.bscf_iter89 = physics_module.bootstrap_fraction_iter89(
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

        betat = (
            physics_variables.beta
            * physics_variables.btot**2
            / physics_variables.bt**2
        )
        current_drive_variables.bscf_nevins = (
            current_drive_variables.cboot
            * physics_module.bootstrap_fraction_nevins(
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
            * physics_module.bootstrap_fraction_wilson(
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
        current_drive_variables.diacf_hender = (
            physics_module.diamagnetic_fraction_hender(physics_variables.beta)
        )

        # SCENE scaling for diamagnetic current
        current_drive_variables.diacf_scene = physics_module.diamagnetic_fraction_scene(
            physics_variables.beta, physics_variables.q95, physics_variables.q0
        )

        # Pfirsch-Schlüter scaling for diamagnetic current
        current_drive_variables.pscf_scene = physics_module.ps_fraction_scene(
            physics_variables.beta
        )

        current_drive_variables.bscf_sauter = (
            current_drive_variables.cboot * physics_module.bootstrap_fraction_sauter()
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
            self.current_drive.cudriv(False)

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
