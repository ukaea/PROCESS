import logging
import numpy
from scipy import integrate
from dataclasses import dataclass
from process.fortran import physics_variables, physics_module, constants
import process.impurity_radiation as impurity


logger = logging.getLogger(__name__)

ATOMIC_MASS_DEUTERIUM = 2.0
ATOMIC_MASS_TRITIUM = 3.0

REACTION_CONSTANTS_DT = dict(
    bg=34.3827,
    mrc2=1.124656e6,
    cc1=1.17302e-9,
    cc2=1.51361e-2,
    cc3=7.51886e-2,
    cc4=4.60643e-3,
    cc5=1.35000e-2,
    cc6=-1.06750e-4,
    cc7=1.36600e-5,
)

REACTION_CONSTANTS_DHE3 = dict(
    bg=68.7508,
    mrc2=1.124572e6,
    cc1=5.51036e-10,
    cc2=6.41918e-3,
    cc3=-2.02896e-3,
    cc4=-1.91080e-5,
    cc5=1.35776e-4,
    cc6=0.0,
    cc7=0.0,
)

REACTION_CONSTANTS_DD1 = dict(
    bg=31.3970,
    mrc2=0.937814e6,
    cc1=5.43360e-12,
    cc2=5.85778e-3,
    cc3=7.68222e-3,
    cc4=0.0,
    cc5=-2.96400e-6,
    cc6=0.0,
    cc7=0.0,
)

REACTION_CONSTANTS_DD2 = dict(
    bg=31.3970,
    mrc2=0.937814e6,
    cc1=5.65718e-12,
    cc2=3.41267e-3,
    cc3=1.99167e-3,
    cc4=0.0,
    cc5=1.05060e-5,
    cc6=0.0,
    cc7=0.0,
)


class FusionReactionRate:
    """Calculate the fusion reaction rate for each reaction case (DT,DHE3,DD1,DD2).
    This subroutine numerically integrates over plasma cross-section to
    find the core plasma fusion power.
    References:
        T & M/PKNIGHT/LOGBOOK24, p.6
    Author:
        P J Knight, CCFE, Culham Science Centre
        G Turkington (UKAEA)
    """

    def __init__(self, plasma_profile):
        self.plasma_profile = plasma_profile
        self.sigvdt = 0.0
        self.pdhe3pv = 0.0
        self.pddpv = 0.0
        self.pdtpv = 0.0
        self.palppv = 0.0
        self.pchargepv = 0.0
        self.pneutpv = 0.0
        self.fusionrate = 0.0
        self.alpharate = 0.0
        self.protonrate = 0.0

    def deuterium_branching(self, ion_temperature: float) -> float:
        """
        Calculate the relative rate of tritium producing D-D reactions to 3He ones based on the volume averaged ion temperature

        Parameters:
            ion_temperature (float): Volume averaged ion temperature in keV

         Returns:
            float: The rate of tritium producing D-D reactions to 3He ones

        Notes:
            - For ion temperatures between 0.5 keV and 200 keV.
            - The deviation of the fit from the R-matrix branching ratio is always smaller than 0.5%.

        References:
            - H.-S. Bosch and G. M. Hale, “Improved formulas for fusion cross-sections and thermal reactivities,”
              Nuclear Fusion, vol. 32, no. 4, pp. 611–631, Apr. 1992,
              doi: https://doi.org/10.1088/0029-5515/32/4/i07.
        ‌
        """
        return (
            1.02934
            - 8.3264e-3 * ion_temperature
            + 1.7631e-4 * ion_temperature**2
            - 1.8201e-6 * ion_temperature**3
            + 6.9855e-9 * ion_temperature**4
        )

    def dt(self):
        """D + T --> 4He + n reaction"""
        dt = BoschHaleConstants(**REACTION_CONSTANTS_DT)
        sigmav = integrate.simpson(
            fint(self.plasma_profile, dt),
            x=self.plasma_profile.neprofile.profile_x,
            dx=self.plasma_profile.neprofile.profile_dx,
        )
        self.sigvdt = sigmav
        etot = 17.59 * constants.electron_charge  # MJ
        # Fusion power produced [MW] per m^3 of plasma
        fusion_power_density = (
            sigmav
            * etot
            * physics_variables.fdeut
            * physics_variables.ftrit
            * physics_variables.deni
            * physics_variables.deni
        )  # MW/m3
        pa = 0.2 * fusion_power_density
        pc = 0.0
        pn = 0.8 * fusion_power_density
        frate = fusion_power_density / etot  # reactions/m3/second
        arate = frate
        prate = 0.0
        self.pdtpv = fusion_power_density
        self.sum_fusion_rates(pa, pc, pn, frate, arate, prate)

    def dhe3(self):
        """D + 3He --> 4He + p reaction"""
        dhe3 = BoschHaleConstants(**REACTION_CONSTANTS_DHE3)
        sigmav = integrate.simpson(
            fint(self.plasma_profile, dhe3),
            x=self.plasma_profile.neprofile.profile_x,
            dx=self.plasma_profile.neprofile.profile_dx,
        )
        etot = 18.35 * constants.electron_charge  # MJ
        # Fusion power produced [MW] per m^3 of plasma
        fusion_power_density = (
            sigmav
            * etot
            * physics_variables.fdeut
            * physics_variables.fhe3
            * physics_variables.deni
            * physics_variables.deni
        )  # MW/m3
        pa = 0.2 * fusion_power_density
        pc = 0.8 * fusion_power_density
        pn = 0.0
        frate = fusion_power_density / etot  # reactions/m3/second
        arate = frate
        prate = frate  # proton production /m3/second
        self.pdhe3pv = fusion_power_density
        self.sum_fusion_rates(pa, pc, pn, frate, arate, prate)

    def dd1(self):
        """D + D --> 3He + n reaction"""
        dd1 = BoschHaleConstants(**REACTION_CONSTANTS_DD1)
        sigmav = integrate.simpson(
            fint(self.plasma_profile, dd1),
            x=self.plasma_profile.neprofile.profile_x,
            dx=self.plasma_profile.neprofile.profile_dx,
        )
        etot = 3.27 * constants.electron_charge  # MJ
        # Fusion power produced [MW] per m^3 of plasma
        fusion_power_density = (
            sigmav
            * etot
            * 0.5
            * physics_variables.fdeut
            * physics_variables.fdeut
            * physics_variables.deni
            * physics_variables.deni
        )  # MW/m3
        pa = 0.0
        pc = 0.25 * fusion_power_density
        pn = 0.75 * fusion_power_density
        frate = fusion_power_density / etot  # reactions/m3/second
        arate = 0.0
        prate = 0.0
        self.pddpv = self.pddpv + fusion_power_density
        self.sum_fusion_rates(pa, pc, pn, frate, arate, prate)

    def dd2(self):
        """D + D --> T + p reaction"""
        dd2 = BoschHaleConstants(**REACTION_CONSTANTS_DD2)
        sigmav = integrate.simpson(
            fint(self.plasma_profile, dd2),
            x=self.plasma_profile.neprofile.profile_x,
            dx=self.plasma_profile.neprofile.profile_dx,
        )
        etot = 4.03 * constants.electron_charge  # MJ
        # Fusion power produced [MW] per m^3 of plasma
        fusion_power_density = (
            sigmav
            * etot
            * 0.5
            * physics_variables.fdeut
            * physics_variables.fdeut
            * physics_variables.deni
            * physics_variables.deni
        )  # MW/m3
        pa = 0.0
        pc = fusion_power_density
        pn = 0.0
        frate = fusion_power_density / etot  # reactions/m3/second
        arate = 0.0
        prate = frate  # proton production /m3/second
        self.pddpv = self.pddpv + fusion_power_density
        self.sum_fusion_rates(pa, pc, pn, frate, arate, prate)

    def sum_fusion_rates(self, pa, pc, pn, frate, arate, prate):
        """Sum the fusion rate at the end of each reaction.

        :param pa: palppv  alpha particle production rate (/m3/s)
        :type pa: float
        :param pc: pchargepv other charged particle fusion power/volume (MW/m3)
        :type pc: float
        :param pn: pneutpv neutron fusion power per volume (MW/m3)
        :type pn: float
        :param frate: fusionrate fusion reaction rate (reactions/m3/s)
        :type frate: float
        :param arate:  alpharate alpha particle fusion power per volume (MW/m3)
        :type arate: float
        :param prate: protonrate proton production rate (/m3/s)
        :type prate: float
        """
        self.palppv = self.palppv + pa
        self.pchargepv = self.pchargepv + pc
        self.pneutpv = self.pneutpv + pn
        self.fusionrate = self.fusionrate + frate
        self.alpharate = self.alpharate + arate
        self.protonrate = self.protonrate + prate

    def calculate_fusion_rates(self):
        """Initiate all the fusion rate calculations."""
        self.dt()
        self.dhe3()
        self.dd1()
        self.dd2()

    def set_physics_variables(self):
        """Set the required physics variables."""
        physics_variables.palppv = self.palppv
        physics_variables.pchargepv = self.pchargepv
        physics_variables.pneutpv = self.pneutpv
        physics_variables.fusionrate = self.fusionrate
        physics_variables.alpharate = self.alpharate
        physics_variables.protonrate = self.protonrate
        physics_module.sigvdt = self.sigvdt
        physics_module.pdtpv = self.pdtpv
        physics_module.pdhe3pv = self.pdhe3pv
        physics_module.pddpv = self.pddpv


def radpwr(plasma_profile):
    """This routine finds the radiation powers in MW/m3 by calling
    relevant routines.
    Author:
        P J Knight, CCFE, Culham Science Centre

    :param plasma_profile:
    :type plasma_profile: PlasmaProfile
    :return: RadpwrData - dataclass returning radiation power
    :rtype: DataClass
    """
    imp_rad = impurity.ImpurityRadiation(plasma_profile)
    imp_rad.calculate_imprad()

    pedgeradpv = imp_rad.radtot - imp_rad.radcore

    # Synchrotron radiation power/volume; assumed to be from core only.
    psyncpv = psync_albajar_fidone()

    # Total core radiation power/volume.
    pcoreradpv = imp_rad.radcore + psyncpv

    # Total radiation power/volume.
    pradpv = imp_rad.radtot + psyncpv

    return RadpwrData(psyncpv, pcoreradpv, pedgeradpv, pradpv)


def psync_albajar_fidone():
    """This routine finds the synchrotron radiation power in MW/m3,
    using the method of Albajar and Fidone.
    References:
        Albajar, Nuclear Fusion 41 (2001) 665
        Fidone, Giruzzi, Granata, Nuclear Fusion 41 (2001) 1755
    Authors:
        P J Knight, CCFE, Culham Science Centre
        R Kemp, CCFE, Culham Science Centre

    :return: psyncpv synchrotron radiation power/volume (MW/m3)
    :rtype: float
    """
    tbet = 2.0e0

    # rpow is the (1-Rsyn) power dependence based on plasma shape
    # (see Fidone)

    rpow = 0.62e0

    kap = 0.0
    de2o = 0.0
    pao = 0.0
    gfun = 0.0
    kfun = 0.0
    dum = 0.0
    psync = 0.0
    psyncpv = 0.0

    kap = physics_variables.vol / (
        2.0e0 * numpy.pi**2 * physics_variables.rmajor * physics_variables.rminor**2
    )

    # No account is taken of pedestal profiles here, other than use of
    # the correct physics_variables.ne0 and physics_variables.te0...

    de2o = 1.0e-20 * physics_variables.ne0
    pao = 6.04e3 * (physics_variables.rminor * de2o) / physics_variables.bt
    gfun = 0.93e0 * (
        1.0e0
        + 0.85e0
        * numpy.exp(-0.82e0 * physics_variables.rmajor / physics_variables.rminor)
    )
    kfun = (physics_variables.alphan + 3.87e0 * physics_variables.alphat + 1.46e0) ** (
        -0.79e0
    )
    kfun = kfun * (1.98e0 + physics_variables.alphat) ** 1.36e0 * tbet**2.14e0
    kfun = kfun * (tbet**1.53e0 + 1.87e0 * physics_variables.alphat - 0.16e0) ** (
        -1.33e0
    )
    dum = (
        1.0e0
        + 0.12e0
        * (physics_variables.te0 / (pao**0.41e0))
        * (1.0e0 - physics_variables.ssync) ** 0.41e0
    )

    # Very high T modification, from Fidone

    dum = dum ** (-1.51e0)

    psync = (
        3.84e-8
        * (1.0e0 - physics_variables.ssync) ** rpow
        * physics_variables.rmajor
        * physics_variables.rminor**1.38e0
    )
    psync = psync * kap**0.79e0 * physics_variables.bt**2.62e0 * de2o**0.38e0
    psync = (
        psync
        * physics_variables.te0
        * (16.0e0 + physics_variables.te0) ** 2.61e0
        * dum
        * gfun
        * kfun
    )

    # psyncpv should be per unit volume; Albajar gives it as total

    psyncpv = psync / physics_variables.vol

    return psyncpv


def fint(plasma_profile, reactionconstants):
    """This function evaluates the integrand for the fusion power
    integration.
    Authors:
        P J Knight, CCFE, Culham Science Centre
        G Turkington (UKAEA)
    References:

    :param plasma_profile: Parameterised temperature and density profiles
    :type plasma_profile: PlasmaProfile
    :param reactionconstants: BoschHale reaction constants
    :type reactionconstants: DataClass
    :return: fint Integrand for the fusion power
    :rtype: numpy.array
    """

    # Since the electron temperature profile is only calculated directly, we scale the ion temperature
    # profile by the ratio of the volume averaged ion to electron temperature
    ion_temperature_profile = (
        physics_variables.ti / physics_variables.te * plasma_profile.teprofile.profile_y
    )

    sigv = numpy.zeros(
        plasma_profile.teprofile.profile_size,
    )

    # Number of fusion reactions per unit volume per particle volume density (m3/s)
    sigv = bosch_hale(ion_temperature_profile, reactionconstants)

    # Integrand for the volume averaged fusion reaction rate sigmav:
    # sigmav = integral(2 rho (sigv(rho) ni(rho)^2) drho),
    # divided by the square of the volume-averaged ion density
    # to retain the dimensions m3/s (thbosch_haleis is multiplied back in later)

    nprof = 1.0 / physics_variables.dene * plasma_profile.neprofile.profile_y
    nprofsq = nprof * nprof
    fint = 2.0 * plasma_profile.teprofile.profile_x * sigv * nprofsq

    return fint


def bosch_hale(temperature_profile, reaction_constants):
    """This routine calculates the volumetric fusion reaction rate
    sigmavgt (m3/s) for one of four nuclear reactions using
    the Bosch-Hale parametrization.
    The valid range of the fit is 0.2 keV < t < 100 keV except for D-3He where it is 0.5 keV < t < 190 keV.
    1 : D-T reaction
    2 : D-3He reaction
    3 : D-D 1st reaction (50% probability)
    4 : D-D 2nd reaction (50% probability)
    Authors:
        R Kemp, CCFE, Culham Science Centre
        P J Knight, CCFE, Culham Science Centre
    References:
        Bosch and Hale, Nuclear Fusion 32 (1992) 611-631

    :param temperature_profile: Plasma temperature profile
    :type temperature_profile: numpy.array
    :param reactionconstants: BoschHale reaction constants
    :type reactionconstants: BoschHaleConstants
    :return: sigmav Volumetric fusion reaction rate
    :rtype: (numpy.array)
    """
    theta1 = (
        temperature_profile
        * (
            reaction_constants.cc2
            + temperature_profile
            * (reaction_constants.cc4 + temperature_profile * reaction_constants.cc6)
        )
        / (
            1.0
            + temperature_profile
            * (
                reaction_constants.cc3
                + temperature_profile
                * (
                    reaction_constants.cc5
                    + temperature_profile * reaction_constants.cc7
                )
            )
        )
    )
    theta = temperature_profile / (1.0 - theta1)

    xi = ((reaction_constants.bg**2) / (4.0 * theta)) ** (1 / 3)

    # Volumetric reaction rate <sigma v> (m3/s)
    # Original form is in [cm3/s], so multiply by 1.0e-6 to convert to [m3/s]
    sigmav = (
        1.0e-6
        * reaction_constants.cc1
        * theta
        * numpy.sqrt(xi / (reaction_constants.mrc2 * temperature_profile**3))
        * numpy.exp(-3.0 * xi)
    )
    # Bosch-Hale also gives value sof maximum deviation of the fit from the input data
    # D-T reaction: max sigmav uncertainty = 0.25%
    # D-3He reaction: max sigmav uncertainty = 2.5%
    # D-D 1st reaction (50% probability)
    # D-D 2nd reaction (50% probability)

    # if t = 0, sigmav = 0. Use this mask to set sigmav to zero.
    t_mask = temperature_profile == 0.0
    sigmav[t_mask] = 0.0

    return sigmav


@dataclass
class BoschHaleConstants:
    """DataClass which holds the constants required for the Bosch Hale calculation
    for a given fusion reaction.
    """

    bg: float
    mrc2: float
    cc1: float
    cc2: float
    cc3: float
    cc4: float
    cc5: float
    cc6: float
    cc7: float


@dataclass
class RadpwrData:
    """DataClass which holds the output of the function radpwr"""

    psyncpv: float
    pcoreradpv: float
    pedgeradpv: float
    pradpv: float


def palph2(
    bp,
    bt,
    dene,
    deni,
    dnitot,
    falpe,
    falpi,
    palpnb,
    pchargepv,
    pneutpv,
    ten,
    tin,
    vol,
    palppv,
    ifalphap,
):
    """
    Fusion power and fast alpha pressure calculations.
    ITER Physics Design Guidelines: 1989 [IPDG89], N. A. Uckan et al,
    ITER Documentation Series No.10, IAEA/ITER/DS/10, IAEA, Vienna, 1990
    D J Ward, UKAEA Fusion: F/PL/PJK/PROCESS/CODE/050

    :param bp: poloidal field (T)
    :param bt: totoidal field on axis (T)
    :param dene: electron density (m^-3)
    :param deni: fuel ion density (m^-3)
    :param dnitot: total ion density (m^-3)
    :param falpe: fraction of alpha energy to electrons
    :param falpi: fraction of alpha energy to ions
    :param palpnb: alpha power from hot neutral beam ions (MW)
    :param pchargepv: other charged particle fusion power/volume (MW/m3)
    :param pneutpv: neutron fusion power per volume (MW/m3)
    :param ten: density-weighted electron temperature (keV)
    :param tin: density-weighted ion temperature (keV)
    :param vol: plasma volume (m3)
    :param palppv: alpha power per volume (MW/m3)
    :param ifalphap: switch for fast alpha pressure method

    :return: neutron fusion power per volume (MW/m3), alpha power (MW),
    neutron fusion power (MW), other charged particle fusion power (MW),
    fast alpha beta component, alpha power per volume (MW/m3),
    alpha power per volume to electrons (MW/m3), alpha power per volume to ions (MW/m3),
    charged particle fusion power (MW), fusion power (MW)
    """

    # Add neutral beam alpha power / volume
    palppv_out = palppv + palpnb / vol

    # Add extra neutron power
    pneutpv_out = pneutpv + 4.0 * palpnb / vol

    # Total alpha power
    palpmw = palppv_out * vol

    # Total non-alpha charged particle power
    pchargemw = pchargepv * vol

    # Total neutron power
    pneutmw = pneutpv_out * vol

    # Total fusion power
    powfmw = palpmw + pneutmw + pchargemw

    # Charged particle fusion power
    pfuscmw = palpmw + pchargemw

    # Alpha power to electrons and ions (used with electron
    # and ion power balance equations only)
    # No consideration of pchargepv here...
    palpipv = physics_variables.falpha * palppv_out * falpi
    palpepv = physics_variables.falpha * palppv_out * falpe

    # Determine average fast alpha density
    if physics_variables.fdeut < 1.0:

        betath = (
            2.0e3
            * constants.rmu0
            * constants.electron_charge
            * (dene * ten + dnitot * tin)
            / (bt**2 + bp**2)
        )

        # jlion: This "fact" model is heavily flawed for smaller temperatures! It is unphysical for a stellarator (high n low T)
        # IPDG89 fast alpha scaling
        if ifalphap == 0:
            fact = min(0.3, 0.29 * (deni / dene) ** 2 * ((ten + tin) / 20.0 - 0.37))

        # Modified scaling, D J Ward
        else:
            fact = min(
                0.30,
                0.26
                * (deni / dene) ** 2
                * numpy.sqrt(max(0.0, ((ten + tin) / 20.0 - 0.65))),
            )

        fact = max(fact, 0.0)
        fact2 = palppv_out / palppv
        betaft = betath * fact * fact2

    else:  # negligible alpha production, palppv = palpnb = 0
        betaft = 0.0

    return (
        pneutpv_out,
        palpmw,
        pneutmw,
        pchargemw,
        betaft,
        palppv_out,
        palpepv,
        palpipv,
        pfuscmw,
        powfmw,
    )


def beamfus(
    beamfus0,
    betbm0,
    bp,
    bt,
    cnbeam,
    dene,
    deni,
    dlamie,
    ealphadt,
    enbeam,
    fdeut,
    ftrit,
    ftritbm,
    sigvdt,
    ten,
    tin,
    vol,
    zeffai,
):
    """Routine to calculate beam slowing down properties
    author: P J Knight, CCFE, Culham Science Centre

    :param beamfus0: multiplier for beam-background fusion calculation
    :param betbm0: leading coefficient for neutral beam beta fraction
    :param bp: poloidal field (T)
    :param bt: toroidal field on axis (T)
    :param cnbeam: neutral beam current (A)
    :param dene: electron density (m^-3)
    :param deni: fuel ion density (m^-3)
    :param dlamie: ion-electron coulomb logarithm
    :param ealphadt: alpha particle birth energy (D-T) (keV)
    :param enbeam: neutral beam energy (keV)
    :param fdeut: deuterium fraction of main plasma
    :param ftrit: tritium fraction of main plasma
    :param ftritbm: tritium fraction of neutral beam
    :param sigvdt: profile averaged <sigma v> for D-T (m3/s)
    :param ten: density-weighted electron temperature (keV)
    :param tin: density-weighted ion temperature (keV)
    :param vol: plasma volume (m3)
    :param zeffai: mass weighted plasma effective charge

    :returns: neutral beam beta component, hot beam ion density (m^-3),
    alpha power from hot neutral beam ions (MW)
    """

    tausl = (
        1.99e19 * (2.0 * (1.0 - ftritbm) + (3.0 * ftritbm)) * (ten**1.5 / dene) / dlamie
    )

    # Critical energy for electron/ion slowing down of the beam ion
    # (deuterium and tritium neutral beams, respectively) (keV)

    ecritd = 14.8 * ten * 2.0 * zeffai**0.6666 * (dlamie + 4.0) / dlamie
    ecritt = ecritd * 1.5

    # Deuterium and tritium ion densities

    denid = deni * fdeut
    denit = deni * ftrit

    palpdb, palptb, dnbeam2, ehotnb = beamcalc(
        denid,
        denit,
        ealphadt,
        enbeam,
        ecritd,
        ecritt,
        tausl,
        ftritbm,
        cnbeam,
        tin,
        vol,
        sigvdt,
    )

    # Neutral beam alpha power

    palpnb = beamfus0 * (palpdb + palptb)

    # Neutral beam beta

    betanb = betbm0 * 4.03e-22 * 0.66666 * dnbeam2 * ehotnb / (bt**2 + bp**2)

    return betanb, dnbeam2, palpnb


def beamcalc(
    nd, nt, ealphadt, ebeam, ecritd, ecritt, tausbme, ftritbm, ibeam, ti, vol, svdt
):
    """Neutral beam alpha power and ion energy
    author: P J Knight, CCFE, Culham Science Centre

    :param nd: thermal deuterium density (/m3)
    :param nt: thermal tritium density   (/m3)
    :param ealphadt: alpha particle birth energy (D-T) (keV)
    :param ebeam: beam energy (keV)
    :param ecritd: critical energy for electron/ion slowing down of
    the beam ion (deuterium neutral beam) (keV)
    :param ecritt: critical energy for beam slowing down (tritium neutral beam) (keV)
    :param ftritbm: beam tritium fraction (0.0 = deuterium beam)
    :param ibeam: beam current (A)
    :param svdt: profile averaged <sigma v> for D-T (m3/s)
    :param tausbme: beam ion slowing down time on electrons (s)
    :param ti: thermal ion temperature (keV)
    :param vol: plasma volume (m3) (95% flux surface)

    :returns: alpha power from deut. beam-background fusion (MW),
    alpha power from trit. beam-background fusion (MW), hot beam ion density (m^-3),
    average hot beam ion energy (keV)
    """

    # D and T beam current fractions
    ifbmd = ibeam * (1.0 - ftritbm)
    ifbmt = ibeam * ftritbm

    ebmratd = ebeam / ecritd
    vcritd = numpy.sqrt(
        2.0
        * constants.electron_charge
        * 1000.0
        * ecritd
        / (constants.proton_mass * ATOMIC_MASS_DEUTERIUM)
    )
    tauseffd = tausbme / 3.0 * numpy.log(1.0 + (ebmratd) ** 1.5)
    nhotmsd = (1.0 - ftritbm) * ibeam * tauseffd / (constants.electron_charge * vol)

    ebmratt = ebeam / ecritt
    vcritt = numpy.sqrt(
        2.0
        * constants.electron_charge
        * 1000.0
        * ecritt
        / (constants.proton_mass * ATOMIC_MASS_TRITIUM)
    )
    tausefft = tausbme / 3.0 * numpy.log(1.0 + (ebmratt) ** 1.5)
    nhotmst = ftritbm * ibeam * tausefft / (constants.electron_charge * vol)

    nhot = nhotmsd + nhotmst
    ndhot = nhotmsd
    nthot = nhotmst

    # Average hot ion energy from Deng & Emmert, UWFDM-718, Jan 87
    vcds = 2.0 * ecritd * constants.electron_charge * 1000.0 / (2.0 * constants.proton_mass)
    vcts = 2.0 * ecritt * constants.electron_charge * 1000.0 / (3.0 * constants.proton_mass)

    s0d = ifbmd / (constants.electron_charge * vol)
    s0t = ifbmt / (constants.electron_charge * vol)

    xcoefd = (
        ATOMIC_MASS_DEUTERIUM
        * constants.proton_mass
        * tausbme
        * vcds
        * s0d
        / (constants.electron_charge * 1000.0 * 3.0)
    )
    xcoeft = (
        ATOMIC_MASS_TRITIUM
        * constants.proton_mass
        * tausbme
        * vcts
        * s0t
        / (constants.electron_charge * 1000.0 * 3.0)
    )

    presd = xcoefd * xbrak(ebeam, ecritd)
    prest = xcoeft * xbrak(ebeam, ecritt)

    ehotd = 1.5 * presd / ndhot
    ehott = 1.5 * prest / nthot
    ehot = (ndhot * ehotd + nthot * ehott) / nhot

    iabm = 2
    svdhotn = 1e-4 * sgvhot(iabm, vcritd, ebeam)
    iabm = 3
    svthotn = 1e-4 * sgvhot(iabm, vcritt, ebeam)

    palfdb = palphabm(ealphadt, ndhot, nt, svdhotn, vol, ti, svdt)
    palftb = palphabm(ealphadt, nthot, nd, svthotn, vol, ti, svdt)

    return palfdb, palftb, nhot, ehot


def xbrak(e0, ec):
    """Hot ion energy parameter
    author: P J Knight, CCFE, Culham Science Centre

    :param e0: neutral beam energy (keV)
    :param ec: critical energy for electron/ion slowing down of the beam ion (keV)
    """
    xcs = e0 / ec
    xc = numpy.sqrt(xcs)

    t1 = xcs / 2.0
    t2 = (numpy.log((xcs + 2.0 * xc + 1.0) / (xcs - xc + 1.0))) / 6.0

    xarg = (2.0 * xc - 1.0) / numpy.sqrt(3.0)
    t3 = (numpy.arctan(xarg)) / numpy.sqrt(3.0)
    t4 = 0.3022999

    return t1 + t2 - t3 - t4


def palphabm(ealphadt, nbm, nblk, sigv, vol, ti, svdt):
    """Alpha power from beam-background fusion
    author: P J Knight, CCFE, Culham Science Centre

    :param ealphadt: alpha particle birth energy (D-T) (keV)
    :param nblk: thermal ion density (/m3)
    :param nbm: hot beam ion density (/m3)
    :param sigv: hot beam fusion reaction rate (m3/s)
    :param vol: plasma volume (m3)
    :param ti: thermal ion temperature (keV)
    :param svdt: profile averaged <sigma v> for D-T (m3/s)
    """

    # [ti] because bosch_hale expects temperature profile
    # so we pass it a profile of length 1
    ratio = svdt / bosch_hale(
        numpy.array([ti]), BoschHaleConstants(**REACTION_CONSTANTS_DT)
    )
    return (
        constants.electron_charge / 1000.0 * nbm * nblk * sigv * ealphadt * vol * ratio.item()
    )


def sgvhot(rmass_ion, vcrx, ebeam):
    """Hot beam fusion reaction rate
    author: P J Knight, CCFE, Culham Science Centre

    :param rmass_ion: relative atomic mass of the ion (of D or T)
    :param vcrx: critical velocity for electron/ion slowing down of the beam ion (m/s)
    :param ebeam: neutral beam energy (keV)
    """
    # Beam velocity

    vbeams = ebeam * constants.electron_charge * 1000.0 * 2.0 / (rmass_ion * constants.proton_mass)
    vbeam = numpy.sqrt(vbeams)

    xv = vbeam / vcrx
    t1 = 3.0 * vcrx / numpy.log(1.0 + (xv**3))

    svint = integrate.quad(
        _hot_beam_fusion_reaction_rate_integrand, 0.0, xv, args=(vcrx,)
    )[0]

    return t1 * svint


def _hot_beam_fusion_reaction_rate_integrand(u, vcritx):
    """Integrand function for the hot beam fusion reaction rate
    author: P J Knight, CCFE, Culham Science Centre

    :param u: ratio of beam velocity to the critical velocity
    """
    t1 = (u**3) / (1.0 + u**3)

    # vcritx : critical velocity for electron/ion slowing down of beam ion (m/s)
    xvc = vcritx * u
    xvcs = xvc * xvc * constants.proton_mass / (constants.electron_charge * 1000.0)
    t2 = _sigbmfus(xvcs)

    return t1 * t2


def _sigbmfus(vrelsq):
    """Fusion reaction cross-section
    author: P J Knight, CCFE, Culham Science Centre

    The functional form of the cross-section is in terms of the equivalent
    deuterium energy, i.e. for a tritium beam at 500 keV the energy
    used in the cross-section function is 333 keV.

    :param vrelsq: square of the speed of the beam ion (keV/amu)
    """
    a1 = 45.95
    a2 = 5.02e4
    a3 = 1.368e-2
    a4 = 1.076
    a5 = 4.09e2

    # Beam kinetic energy

    ebm = 0.5 * ATOMIC_MASS_DEUTERIUM * vrelsq

    # Set limits on cross-section at low and high beam energies

    if ebm < 10.0:
        return 1.0e-27
    elif ebm > 1.0e4:
        return 8.0e-26
    else:
        t1 = a2 / (1.0e0 + (a3 * ebm - a4) ** 2) + a5
        t2 = ebm * (numpy.exp(a1 / numpy.sqrt(ebm)) - 1.0)
        return 1.0e-24 * t1 / t2
