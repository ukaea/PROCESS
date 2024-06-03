import logging
import numpy
from scipy import integrate
from process.fortran import constants
from dataclasses import dataclass
from process.fortran import physics_variables
from process.fortran import physics_module
import process.impurity_radiation as impurity


logger = logging.getLogger(__name__)


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

    def dt(self):
        """D + T --> 4He + n reaction"""
        dt = BoschHaleConstants(
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
        sigmav = integrate.simpson(
            fint(self.plasma_profile, dt),
            x=self.plasma_profile.neprofile.profile_x,
            dx=self.plasma_profile.neprofile.profile_dx,
        )
        self.sigvdt = sigmav
        etot = 17.59 * constants.echarge  # MJ
        fpow = (
            1.0
            * sigmav
            * etot
            * physics_variables.fdeut
            * physics_variables.ftrit
            * physics_variables.deni
            * physics_variables.deni
        )  # MW/m3
        pa = 0.2 * fpow
        pc = 0.0
        pn = 0.8 * fpow
        frate = fpow / etot  # reactions/m3/second
        arate = frate
        prate = 0.0
        self.pdtpv = fpow
        self.sum_fusion_rates(pa, pc, pn, frate, arate, prate)

    def dhe3(self):
        """D + 3He --> 4He + p reaction"""
        dhe3 = BoschHaleConstants(
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
        sigmav = integrate.simpson(
            fint(self.plasma_profile, dhe3),
            x=self.plasma_profile.neprofile.profile_x,
            dx=self.plasma_profile.neprofile.profile_dx,
        )
        etot = 18.35 * constants.echarge  # MJ
        fpow = (
            1.0
            * sigmav
            * etot
            * physics_variables.fdeut
            * physics_variables.fhe3
            * physics_variables.deni
            * physics_variables.deni
        )  # MW/m3
        pa = 0.2 * fpow
        pc = 0.8 * fpow
        pn = 0.0
        frate = fpow / etot  # reactions/m3/second
        arate = frate
        prate = frate  # proton production /m3/second
        self.pdhe3pv = fpow
        self.sum_fusion_rates(pa, pc, pn, frate, arate, prate)

    def dd1(self):
        """D + D --> 3He + n reaction"""
        dd1 = BoschHaleConstants(
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
        sigmav = integrate.simpson(
            fint(self.plasma_profile, dd1),
            x=self.plasma_profile.neprofile.profile_x,
            dx=self.plasma_profile.neprofile.profile_dx,
        )
        etot = 3.27 * constants.echarge  # MJ
        fpow = (
            1.0
            * sigmav
            * etot
            * 0.5
            * physics_variables.fdeut
            * physics_variables.fdeut
            * physics_variables.deni
            * physics_variables.deni
        )  # MW/m3
        pa = 0.0
        pc = 0.25 * fpow
        pn = 0.75 * fpow
        frate = fpow / etot  # reactions/m3/second
        arate = 0.0
        prate = 0.0
        self.pddpv = self.pddpv + fpow
        self.sum_fusion_rates(pa, pc, pn, frate, arate, prate)

    def dd2(self):
        """D + D --> T + p reaction"""
        dd2 = BoschHaleConstants(
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
        sigmav = integrate.simpson(
            fint(self.plasma_profile, dd2),
            x=self.plasma_profile.neprofile.profile_x,
            dx=self.plasma_profile.neprofile.profile_dx,
        )
        etot = 4.03 * constants.echarge  # MJ
        fpow = (
            1.0
            * sigmav
            * etot
            * 0.5
            * physics_variables.fdeut
            * physics_variables.fdeut
            * physics_variables.deni
            * physics_variables.deni
        )  # MW/m3
        pa = 0.0
        pc = fpow
        pn = 0.0
        frate = fpow / etot  # reactions/m3/second
        arate = 0.0
        prate = frate  # proton production /m3/second
        self.pddpv = self.pddpv + fpow
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

    #  Synchrotron radiation power/volume; assumed to be from core only.
    psyncpv = psync_albajar_fidone()

    #  Total core radiation power/volume.
    pcoreradpv = imp_rad.radcore + psyncpv

    #  Total radiation power/volume.
    pradpv = imp_rad.radtot + psyncpv

    return RadpwrData(
        psyncpv, pcoreradpv, pedgeradpv, pradpv
    )


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

    #  rpow is the (1-Rsyn) power dependence based on plasma shape
    #  (see Fidone)

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

    #  No account is taken of pedestal profiles here, other than use of
    #  the correct physics_variables.ne0 and physics_variables.te0...

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

    #  Very high T modification, from Fidone

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

    #  psyncpv should be per unit volume; Albajar gives it as total

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

    # Local ion temperature (keV) at r/a = rho
    sigv = numpy.zeros(
        plasma_profile.teprofile.profile_size,
    )

    tiofr = (
        physics_variables.ti / physics_variables.te * plasma_profile.teprofile.profile_y
    )

    # Fusion reaction rate (m3/s)
    sigv = bosch_hale(tiofr, reactionconstants)

    # Integrand for the volume averaged fusion reaction rate sigmav:
    # sigmav = integral(2 rho (sigv(rho) ni(rho)^2) drho),
    # divided by the square of the volume-averaged ion density
    # to retain the dimensions m3/s (thbosch_haleis is multiplied back in later)

    nprof = 1.0 / physics_variables.dene * plasma_profile.neprofile.profile_y
    nprofsq = nprof * nprof
    fint = 2.0 * plasma_profile.teprofile.profile_x * sigv * nprofsq

    return fint


def bosch_hale(t, reactionconstants):
    """This routine calculates the volumetric fusion reaction rate
    sigmavgt (m3/s) for one of four nuclear reactions using
    the Bosch-Hale parametrization.
    The valid range of the fit is 0.2 keV < t < 100 keV
    1 : D-T reaction
    2 : D-3He reaction
    3 : D-D 1st reaction (50% probability)
    4 : D-D 2nd reaction (50% probability)
    Authors:
        R Kemp, CCFE, Culham Science Centre
        P J Knight, CCFE, Culham Science Centre
    References:
        Bosch and Hale, Nuclear Fusion 32 (1992) 611-631

    :param t: Plasma temperature profile
    :type t: numpy.array
    :param reactionconstants: BoschHale reaction constants
    :type reactionconstants: BoschHaleConstants
    :return: sigmav Volumetric fusion reaction rate
    :rtype: (numpy.array)
    """
    theta1 = (
        t
        * (
            reactionconstants.cc2
            + t * (reactionconstants.cc4 + t * reactionconstants.cc6)
        )
        / (
            1.0
            + t
            * (
                reactionconstants.cc3
                + t * (reactionconstants.cc5 + t * reactionconstants.cc7)
            )
        )
    )
    theta = t / (1.0 - theta1)

    xi = ((reactionconstants.bg**2) / (4.0 * theta)) ** 0.3333333333

    # Volumetric reaction rate <sigma v> (m3/s)
    sigmav = (
        1.0e-6
        * reactionconstants.cc1
        * theta
        * numpy.sqrt(xi / (reactionconstants.mrc2 * t**3))
        * numpy.exp(-3.0 * xi)
    )
    # if t = 0, sigmav = 0. Use this mask to set sigmav to zero.
    t_mask = t == 0.0
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
