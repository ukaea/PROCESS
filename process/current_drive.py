import logging

import numpy as np

from process import constants
from process import (
    process_output as po,
)
from process.data_structure import (
    current_drive_variables,
    heat_transport_variables,
    physics_variables,
)
from process.exceptions import ProcessError, ProcessValueError
from process.plasma_profiles import PlasmaProfile

logger = logging.getLogger(__name__)


class NeutralBeam:
    def __init__(self, plasma_profile: PlasmaProfile):
        self.outfile = constants.NOUT
        self.plasma_profile = plasma_profile

    def iternb(self):
        """Routine to calculate ITER Neutral Beam current drive parameters
        author: P J Knight, CCFE, Culham Science Centre
        effnbss : output real : neutral beam current drive efficiency (A/W)
        f_p_beam_injected_ions   : output real : fraction of NB power given to ions
        fshine  : output real : shine-through fraction of beam
        This routine calculates the current drive parameters for a
        neutral beam system, based on the 1990 ITER model.
        ITER Physics Design Guidelines: 1989 [IPDG89], N. A. Uckan et al,
        ITER Documentation Series No.10, IAEA/ITER/DS/10, IAEA, Vienna, 1990
        """
        # Check argument sanity
        if (
            1 + physics_variables.eps
        ) < current_drive_variables.f_radius_beam_tangency_rmajor:
            raise ProcessValueError(
                "Imminent negative square root argument; NBI will miss plasma completely",
                eps=physics_variables.eps,
                f_radius_beam_tangency_rmajor=current_drive_variables.f_radius_beam_tangency_rmajor,
            )

        # Calculate beam path length to centre
        dpath = physics_variables.rmajor * np.sqrt(
            (1.0 + physics_variables.eps) ** 2
            - current_drive_variables.f_radius_beam_tangency_rmajor**2
        )

        # Calculate beam stopping cross-section
        sigstop = self.sigbeam(
            current_drive_variables.e_beam_kev / physics_variables.m_beam_amu,
            physics_variables.temp_plasma_electron_vol_avg_kev,
            physics_variables.nd_plasma_electrons_vol_avg,
            physics_variables.f_nd_alpha_electron,
            physics_variables.f_nd_plasma_carbon_electron,
            physics_variables.f_nd_plasma_oxygen_electron,
            physics_variables.f_nd_plasma_iron_argon_electron,
        )

        # Calculate number of decay lengths to centre
        current_drive_variables.n_beam_decay_lengths_core = (
            dpath * physics_variables.nd_plasma_electrons_vol_avg * sigstop
        )

        # Shine-through fraction of beam
        fshine = np.exp(
            -2.0 * dpath * physics_variables.nd_plasma_electrons_vol_avg * sigstop
        )
        fshine = max(fshine, 1e-20)

        # Deuterium and tritium beam densities
        dend = physics_variables.nd_plasma_fuel_ions_vol_avg * (
            1.0 - current_drive_variables.f_beam_tritium
        )
        dent = (
            physics_variables.nd_plasma_fuel_ions_vol_avg
            * current_drive_variables.f_beam_tritium
        )

        # Power split to ions / electrons
        f_p_beam_injected_ions = self.cfnbi(
            physics_variables.m_beam_amu,
            current_drive_variables.e_beam_kev,
            physics_variables.temp_plasma_electron_density_weighted_kev,
            physics_variables.nd_plasma_electrons_vol_avg,
            dend,
            dent,
            physics_variables.n_charge_plasma_effective_mass_weighted_vol_avg,
            physics_variables.dlamie,
        )

        # Current drive efficiency
        effnbss = current_drive_variables.f_radius_beam_tangency_rmajor * self.etanb(
            physics_variables.m_beam_amu,
            physics_variables.alphan,
            physics_variables.alphat,
            physics_variables.aspect,
            physics_variables.nd_plasma_electrons_vol_avg,
            current_drive_variables.e_beam_kev,
            physics_variables.rmajor,
            physics_variables.temp_plasma_electron_density_weighted_kev,
            physics_variables.n_charge_plasma_effective_vol_avg,
        )

        return effnbss, f_p_beam_injected_ions, fshine

    def culnbi(self):
        """Routine to calculate Neutral Beam current drive parameters
        author: P J Knight, CCFE, Culham Science Centre
        effnbss : output real : neutral beam current drive efficiency (A/W)
        f_p_beam_injected_ions   : output real : fraction of NB power given to ions
        fshine  : output real : shine-through fraction of beam
        This routine calculates Neutral Beam current drive parameters
        using the corrections outlined in AEA FUS 172 to the ITER method.
        <P>The result cannot be guaranteed for devices with aspect ratios far
        from that of ITER (approx. 2.8).
        AEA FUS 172: Physics Assessment for the European Reactor Study
        """
        if (
            1.0e0 + physics_variables.eps
        ) < current_drive_variables.f_radius_beam_tangency_rmajor:
            raise ProcessValueError(
                "Imminent negative square root argument; NBI will miss plasma completely",
                eps=physics_variables.eps,
                f_radius_beam_tangency_rmajor=current_drive_variables.f_radius_beam_tangency_rmajor,
            )

        #  Calculate beam path length to centre

        dpath = physics_variables.rmajor * np.sqrt(
            (1.0e0 + physics_variables.eps) ** 2
            - current_drive_variables.f_radius_beam_tangency_rmajor**2
        )

        #  Calculate beam stopping cross-section

        sigstop = self.sigbeam(
            current_drive_variables.e_beam_kev / physics_variables.m_beam_amu,
            physics_variables.temp_plasma_electron_vol_avg_kev,
            physics_variables.nd_plasma_electrons_vol_avg,
            physics_variables.f_nd_alpha_electron,
            physics_variables.f_nd_plasma_carbon_electron,
            physics_variables.f_nd_plasma_oxygen_electron,
            physics_variables.f_nd_plasma_iron_argon_electron,
        )

        #  Calculate number of decay lengths to centre

        current_drive_variables.n_beam_decay_lengths_core = (
            dpath * physics_variables.nd_plasma_electron_line * sigstop
        )

        #  Shine-through fraction of beam

        fshine = np.exp(
            -2.0e0 * dpath * physics_variables.nd_plasma_electron_line * sigstop
        )
        fshine = max(fshine, 1.0e-20)

        #  Deuterium and tritium beam densities

        dend = physics_variables.nd_plasma_fuel_ions_vol_avg * (
            1.0e0 - current_drive_variables.f_beam_tritium
        )
        dent = (
            physics_variables.nd_plasma_fuel_ions_vol_avg
            * current_drive_variables.f_beam_tritium
        )

        #  Power split to ions / electrons

        f_p_beam_injected_ions = self.cfnbi(
            physics_variables.m_beam_amu,
            current_drive_variables.e_beam_kev,
            physics_variables.temp_plasma_electron_density_weighted_kev,
            physics_variables.nd_plasma_electrons_vol_avg,
            dend,
            dent,
            physics_variables.n_charge_plasma_effective_mass_weighted_vol_avg,
            physics_variables.dlamie,
        )

        #  Current drive efficiency

        effnbss = self.etanb2(
            physics_variables.m_beam_amu,
            physics_variables.alphan,
            physics_variables.alphat,
            physics_variables.aspect,
            physics_variables.nd_plasma_electrons_vol_avg,
            physics_variables.nd_plasma_electron_line,
            current_drive_variables.e_beam_kev,
            current_drive_variables.f_radius_beam_tangency_rmajor,
            fshine,
            physics_variables.rmajor,
            physics_variables.rminor,
            physics_variables.temp_plasma_electron_density_weighted_kev,
            physics_variables.n_charge_plasma_effective_vol_avg,
        )

        return effnbss, f_p_beam_injected_ions, fshine

    def etanb2(
        self,
        m_beam_amu,
        alphan,
        alphat,
        aspect,
        nd_plasma_electrons_vol_avg,
        nd_plasma_electron_line,
        e_beam_kev,
        f_radius_beam_tangency_rmajor,
        fshine,
        rmajor,
        rminor,
        temp_plasma_electron_density_weighted_kev,
        zeff,
    ):
        """Routine to find neutral beam current drive efficiency
        using the ITER 1990 formulation, plus correction terms
        outlined in Culham Report AEA FUS 172
        author: P J Knight, CCFE, Culham Science Centre
        m_beam_amu   : input real : beam ion mass (amu)
        alphan  : input real : density profile factor
        alphat  : input real : temperature profile factor
        aspect  : input real : aspect ratio
        nd_plasma_electrons_vol_avg    : input real : volume averaged electron density (m**-3)
        nd_plasma_electron_line    : input real : line averaged electron density (m**-3)
        e_beam_kev  : input real : neutral beam energy (keV)
        f_radius_beam_tangency_rmajor  : input real : R_tangent / R_major for neutral beam injection
        fshine  : input real : shine-through fraction of beam
        rmajor  : input real : plasma major radius (m)
        rminor  : input real : plasma minor radius (m)
        temp_plasma_electron_density_weighted_kev     : input real : density weighted average electron temperature (keV)
        zeff    : input real : plasma effective charge
        This routine calculates the current drive efficiency in A/W of
        a neutral beam system, based on the 1990 ITER model,
        plus correction terms outlined in Culham Report AEA FUS 172.
        <P>The formulae are from AEA FUS 172, unless denoted by IPDG89.
        AEA FUS 172: Physics Assessment for the European Reactor Study
        ITER Physics Design Guidelines: 1989 [IPDG89], N. A. Uckan et al,
        ITER Documentation Series No.10, IAEA/ITER/DS/10, IAEA, Vienna, 1990
        """
        #  Charge of beam ions
        zbeam = 1.0

        #  Fitting factor (IPDG89)
        bbd = 1.0

        #  Volume averaged electron density (10**20 m**-3)
        dene20 = nd_plasma_electrons_vol_avg / 1e20

        #  Line averaged electron density (10**20 m**-3)
        dnla20 = nd_plasma_electron_line / 1e20

        #  Critical energy (MeV) (power to electrons = power to ions) (IPDG89)
        #  N.B. temp_plasma_electron_density_weighted_kev is in keV
        ecrit = 0.01 * m_beam_amu * temp_plasma_electron_density_weighted_kev

        #  Beam energy in MeV
        ebmev = e_beam_kev / 1e3

        #  x and y coefficients of function J0(x,y) (IPDG89)
        xjs = ebmev / (bbd * ecrit)
        xj = np.sqrt(xjs)

        yj = 0.8 * zeff / m_beam_amu

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
        r = max(rmajor, rmajor * f_radius_beam_tangency_rmajor)
        eps1 = rminor / r

        if (1.0 + eps1) < f_radius_beam_tangency_rmajor:
            raise ProcessValueError(
                "Imminent negative square root argument; NBI will miss plasma completely",
                eps=eps1,
                f_radius_beam_tangency_rmajor=f_radius_beam_tangency_rmajor,
            )

        d = rmajor * np.sqrt((1.0 + eps1) ** 2 - f_radius_beam_tangency_rmajor**2)

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
        gamnb = (
            5.0
            * abd
            * 0.1
            * temp_plasma_electron_density_weighted_kev
            * (1.0 - fshine)
            * f_radius_beam_tangency_rmajor
            * j0
            / 0.2
            * ffac
        )

        #  Current drive efficiency (A/W)
        return gamnb / (dene20 * rmajor)

    def etanb(
        self,
        m_beam_amu,
        alphan,
        alphat,
        aspect,
        nd_plasma_electrons_vol_avg,
        ebeam,
        rmajor,
        temp_plasma_electron_density_weighted_kev,
        zeff,
    ):
        """Routine to find neutral beam current drive efficiency
        using the ITER 1990 formulation
        author: P J Knight, CCFE, Culham Science Centre
        m_beam_amu   : input real : beam ion mass (amu)
        alphan  : input real : density profile factor
        alphat  : input real : temperature profile factor
        aspect  : input real : aspect ratio
        nd_plasma_electrons_vol_avg    : input real : volume averaged electron density (m**-3)
        ebeam  : input real : neutral beam energy (keV)
        rmajor  : input real : plasma major radius (m)
        temp_plasma_electron_density_weighted_kev     : input real : density weighted average electron temp. (keV)
        zeff    : input real : plasma effective charge
        This routine calculates the current drive efficiency of
        a neutral beam system, based on the 1990 ITER model.
        ITER Physics Design Guidelines: 1989 [IPDG89], N. A. Uckan et al,
        ITER Documentation Series No.10, IAEA/ITER/DS/10, IAEA, Vienna, 1990
        """

        zbeam = 1.0
        bbd = 1.0

        dene20 = 1e-20 * nd_plasma_electrons_vol_avg

        # Ratio of E_beam/E_crit
        xjs = ebeam / (
            bbd * 10.0 * m_beam_amu * temp_plasma_electron_density_weighted_kev
        )
        xj = np.sqrt(xjs)

        yj = 0.8 * zeff / m_beam_amu

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

        return (
            abd
            * (5.0 / rmajor)
            * (0.1 * temp_plasma_electron_density_weighted_kev / dene20)
            * rjfunc
            / 0.2
            * ffac
        )

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
        a = np.array([
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
        ])

        b = np.array([
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
        ])

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

    def cfnbi(
        self,
        afast,
        efast,
        te,
        ne,
        _nd,
        _nt,
        n_charge_plasma_effective_mass_weighted_vol_avg,
        xlmbda,
    ):
        """Routine to calculate the fraction of the fast particle energy
        coupled to the ions
        author: P J Knight, CCFE, Culham Science Centre
        afast   : input real : mass of fast particle (units of proton mass)
        efast   : input real : energy of fast particle (keV)
        te      : input real : density weighted average electron temp. (keV)
        ne      : input real : volume averaged electron density (m**-3)
        nd      : input real : deuterium beam density (m**-3)
        nt      : input real : tritium beam density (m**-3)
        n_charge_plasma_effective_mass_weighted_vol_avg  : input real : mass weighted plasma effective charge
        xlmbda  : input real : ion-electron coulomb logarithm
        f_p_beam_injected_ions   : output real : fraction of fast particle energy coupled to ions
        This routine calculates the fast particle energy coupled to
        the ions in the neutral beam system.
        """
        # atmd = 2.0
        atmdt = 2.5
        # atmt = 3.0
        c = 3.0e8
        me = constants.ELECTRON_MASS
        # zd = 1.0
        # zt = 1.0

        # xlbd = self.xlmbdabi(afast, atmd, efast, te, ne)
        # xlbt = self.xlmbdabi(afast, atmt, efast, te, ne)

        # sum = nd * zd * zd * xlbd / atmd + nt * zt * zt * xlbt / atmt
        # ecritfix = 16.0e0 * te * afast * (sum / (ne * xlmbda)) ** (2.0e0 / 3.0e0)

        xlmbdai = self.xlmbdabi(afast, atmdt, efast, te, ne)
        sumln = n_charge_plasma_effective_mass_weighted_vol_avg * xlmbdai / xlmbda
        xlnrat = (
            3.0e0 * np.sqrt(np.pi) / 4.0e0 * me / constants.PROTON_MASS * sumln
        ) ** (2.0e0 / 3.0e0)
        ve = c * np.sqrt(2.0e0 * te / 511.0e0)

        ecritfi = (
            afast
            * constants.PROTON_MASS
            * ve
            * ve
            * xlnrat
            / (2.0e0 * constants.ELECTRON_CHARGE * 1.0e3)
        )

        x = np.sqrt(efast / ecritfi)
        t1 = np.log((x * x - x + 1.0e0) / ((x + 1.0e0) ** 2))
        thx = (2.0e0 * x - 1.0e0) / np.sqrt(3.0e0)
        t2 = 2.0e0 * np.sqrt(3.0e0) * (np.arctan(thx) + np.pi / 6.0e0)

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


class ElectronCyclotron:
    def __init__(self, plasma_profile: PlasmaProfile):
        self.outfile = constants.NOUT
        self.plasma_profile = plasma_profile

    def culecd(self):
        """Routine to calculate Electron Cyclotron current drive efficiency
        author: M R O'Brien, CCFE, Culham Science Centre
        author: P J Knight, CCFE, Culham Science Centre
        effrfss : output real : electron cyclotron current drive efficiency (A/W)
        This routine calculates the current drive parameters for a
        electron cyclotron system, based on the AEA FUS 172 model.
        AEA FUS 172: Physics Assessment for the European Reactor Study
        """
        rrr = 1.0e0 / 3.0e0

        #  Temperature
        tlocal = self.plasma_profile.teprofile.calculate_profile_y(
            rrr,
            physics_variables.radius_plasma_pedestal_temp_norm,
            physics_variables.temp_plasma_electron_on_axis_kev,
            physics_variables.temp_plasma_pedestal_kev,
            physics_variables.temp_plasma_separatrix_kev,
            physics_variables.alphat,
            physics_variables.tbeta,
        )

        #  Density (10**20 m**-3)
        dlocal = 1.0e-20 * self.plasma_profile.neprofile.calculate_profile_y(
            rrr,
            physics_variables.radius_plasma_pedestal_density_norm,
            physics_variables.nd_plasma_electron_on_axis,
            physics_variables.nd_plasma_pedestal_electron,
            physics_variables.nd_plasma_separatrix_electron,
            physics_variables.alphan,
        )

        #  Inverse aspect ratio
        epsloc = rrr * physics_variables.rminor / physics_variables.rmajor

        #  Effective charge (use average value)
        zlocal = physics_variables.n_charge_plasma_effective_vol_avg

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
        AEA FUS 172: Physics Assessment for the European Reactor Study
        ITER Physics Design Guidelines: 1989 [IPDG89], N. A. Uckan et al,
        ITER Documentation Series No.10, IAEA/ITER/DS/10, IAEA, Vienna, 1990
        """
        mcsq = (
            constants.ELECTRON_MASS * 2.9979e8**2 / (1.0e3 * constants.ELECTRON_VOLT)
        )  # keV
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
            raise ProcessValueError("Negative normalised current drive efficiency")
        return ecgam

    def electron_cyclotron_fenstermacher(
        self,
        temp_plasma_electron_density_weighted_kev: float,
        rmajor: float,
        dene20: float,
        dlamee: float,
    ) -> float:
        """
        Routine to calculate Fenstermacher Electron Cyclotron heating efficiency.

        :param temp_plasma_electron_density_weighted_kev: Density weighted average electron temperature keV.
        :type temp_plasma_electron_density_weighted_kev: float
        :param zeff: Plasma effective charge.
        :type zeff: float
        :param rmajor: Major radius of the plasma in meters.
        :type rmajor: float
        :param dene20: Volume averaged electron density in 1x10^20 m^-3.
        :type dene20: float
        :param dlamee: Electron collision frequency in 1/s.
        :type dlamee: float

        :return: The calculated electron cyclotron heating efficiency in A/W.
        :rtype: float

        :notes:

        :references:
            - T.C. Hender et al., 'Physics Assessment of the European Reactor Study', AEA FUS 172, 1992.

        """

        return (0.21e0 * temp_plasma_electron_density_weighted_kev) / (
            rmajor * dene20 * dlamee
        )

    def electron_cyclotron_freethy(
        self,
        te: float,
        zeff: float,
        rmajor: float,
        nd_plasma_electrons_vol_avg: float,
        b_plasma_toroidal_on_axis: float,
        n_ecrh_harmonic: int,
        i_ecrh_wave_mode: int,
    ) -> float:
        """
        Calculate the Electron Cyclotron current drive efficiency using the Freethy model.

        This function computes the ECCD efficiency based on the electron temperature,
        effective charge, major radius, electron density, magnetic field, harmonic number,
        and wave mode.

        :param te: Volume averaged electron temperature in keV.
        :type te: float
        :param zeff: Plasma effective charge.
        :type zeff: float
        :param rmajor: Major radius of the plasma in meters.
        :type rmajor: float
        :param nd_plasma_electrons_vol_avg: Volume averaged electron density in m^-3.
        :type nd_plasma_electrons_vol_avg: float
        :param b_plasma_toroidal_on_axis: Toroidal magnetic field in Tesla.
        :type b_plasma_toroidal_on_axis: float
        :param n_ecrh_harmonic: Cyclotron harmonic number (fundamental used as default).
        :type n_ecrh_harmonic: int
        :param i_ecrh_wave_mode: Wave mode switch (0 for O-mode, 1 for X-mode).
        :type i_ecrh_wave_mode: int

        :return: The calculated absolute ECCD efficiency in A/W.
        :rtype: float

        :notes:
            - Plasma coupling only occurs if the plasma cut-off is below the cyclotron harmonic.
            - The density factor accounts for this behavior.

        :references:
            - Freethy, S., PROCESS issue #2994.
        """

        # Cyclotron frequency
        fc = (
            1
            / (2 * np.pi)
            * constants.ELECTRON_CHARGE
            * b_plasma_toroidal_on_axis
            / constants.ELECTRON_MASS
        )

        # Plasma frequency
        fp = (
            1
            / (2 * np.pi)
            * np.sqrt(
                (nd_plasma_electrons_vol_avg)
                * constants.ELECTRON_CHARGE**2
                / (constants.ELECTRON_MASS * constants.EPSILON0)
            )
        )

        # Scaling factor for ECCD efficiency
        xi_CD = 0.18e0  # Tuned to the results of a GRAY study
        xi_CD *= 4.8e0 / (2 + zeff)  # Zeff correction

        # ECCD efficiency
        eta_cd = xi_CD * te / (3.27e0 * rmajor * (nd_plasma_electrons_vol_avg / 1.0e19))

        # Determine the cut-off frequency based on wave mode
        if i_ecrh_wave_mode == 0:  # O-mode case
            f_cutoff = fp
        elif i_ecrh_wave_mode == 1:  # X-mode case
            f_cutoff = 0.5 * (fc + np.sqrt(n_ecrh_harmonic * fc**2 + 4 * fp**2))
        else:
            raise ValueError("Invalid wave mode. Use 0 for O-mode or 1 for X-mode.")

        # Plasma coupling factor
        a = 0.1  # Controls sharpness of the transition
        cutoff_factor = 0.5 * (
            1 + np.tanh((2 / a) * ((n_ecrh_harmonic * fc - f_cutoff) / fp - a))
        )

        # Final ECCD efficiency
        return eta_cd * cutoff_factor

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
            raise ProcessValueError("Invalid argument", arg=arg)

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
            raise ProcessError("legend: Solution has not converged")


class IonCyclotron:
    def __init__(self, plasma_profile: PlasmaProfile):
        self.outfile = constants.NOUT
        self.plasma_profile = plasma_profile

    def ion_cyclotron_ipdg89(
        self,
        temp_plasma_electron_density_weighted_kev: float,
        zeff: float,
        rmajor: float,
        dene20: float,
    ) -> float:
        """
        Routine to calculate IPDG89 Ion Cyclotron heating efficiency.

        This function computes the ion cyclotron heating efficiency based on
        the electron temperature, effective charge, major radius, and electron density.

        :param temp_plasma_electron_density_weighted_kev: Density weighted average electron temperature keV.
        :type temp_plasma_electron_density_weighted_kev: float
        :param zeff: Plasma effective charge.
        :type zeff: float
        :param rmajor: Major radius of the plasma in meters.
        :type rmajor: float
        :param nd_plasma_electrons_vol_avg: Volume averaged electron density in 1x10^20 m^-3.
        :type nd_plasma_electrons_vol_avg: float

        :return: The calculated ion cyclotron heating efficiency in A/W.
        :rtype: float

        :notes:
        - The 0.1 term is to convert the temperature into 10 keV units
        - The original formula is for the normalised current drive efficnecy
          hence the addition of the density and majro radius terms to get back to an absolute value

        :references:
        - N.A. Uckan and ITER Physics Group, 'ITER Physics Design Guidelines: 1989',
          https://inis.iaea.org/collection/NCLCollectionStore/_Public/21/068/21068960.pdf

        - T.C. Hender et al., 'Physics Assessment of the European Reactor Study', AEA FUS 172, 1992.
        """

        return (
            (0.63e0 * 0.1e0 * temp_plasma_electron_density_weighted_kev)
            / (2.0e0 + zeff)
        ) / (rmajor * dene20)


class ElectronBernstein:
    def __init__(self, plasma_profile: PlasmaProfile):
        self.outfile = constants.NOUT
        self.plasma_profile = plasma_profile

    def electron_bernstein_freethy(
        self,
        te: float,
        rmajor: float,
        dene20: float,
        b_plasma_toroidal_on_axis: float,
        n_ecrh_harmonic: int,
        xi_ebw: float,
    ) -> float:
        """
        Calculate the Electron Bernstein Wave (EBW) current drive efficiency using the Freethy model.

        This function computes the EBW current drive efficiency based on the electron temperature,
        major radius, electron density, magnetic field, harmonic number, and scaling factor.

        :param te: Volume averaged electron temperature in keV.
        :type te: float
        :param rmajor: Major radius of the plasma in meters.
        :type rmajor: float
        :param dene20: Volume averaged electron density in units of 10^20 m^-3.
        :type dene20: float
        :param b_plasma_toroidal_on_axis: Toroidal magnetic field in Tesla.
        :type b_plasma_toroidal_on_axis: float
        :param n_ecrh_harmonic: Cyclotron harmonic number (fundamental used as default).
        :type n_ecrh_harmonic: int
        :param xi_ebw: Scaling factor for EBW efficiency.
        :type xi_ebw: float

        :return: The calculated absolute EBW current drive efficiency in A/W.
        :rtype: float

        :notes:
            - EBWs can only couple to plasma if the cyclotron harmonic is above the plasma density cut-off.
            - The density factor accounts for this behavior.

        :references:
            - Freethy, S., PROCESS issue #1262.
        """

        # Normalised current drive efficiency gamma
        eta_cd_norm = (xi_ebw / 32.7e0) * te

        # Absolute current drive efficiency
        eta_cd = eta_cd_norm / (dene20 * rmajor)

        # EBWs can only couple to plasma if cyclotron harmonic is above plasma density cut-off;
        # this behavior is captured in the following function:
        # constant 'a' controls sharpness of transition
        a = 0.1e0

        fc = (
            1.0e0
            / (2.0e0 * np.pi)
            * n_ecrh_harmonic
            * constants.ELECTRON_CHARGE
            * b_plasma_toroidal_on_axis
            / constants.ELECTRON_MASS
        )

        fp = (
            1.0e0
            / (2.0e0 * np.pi)
            * np.sqrt(
                dene20
                * 1.0e20
                * constants.ELECTRON_CHARGE**2
                / (constants.ELECTRON_MASS * constants.EPSILON0)
            )
        )

        density_factor = 0.5e0 * (1.0e0 + np.tanh((2.0e0 / a) * ((fp - fc) / fp - a)))

        return eta_cd * density_factor


class LowerHybrid:
    def __init__(self, plasma_profile: PlasmaProfile):
        self.outfile = constants.NOUT
        self.plasma_profile = plasma_profile

    def cullhy(self):
        """Routine to calculate Lower Hybrid current drive efficiency
        author: P J Knight, CCFE, Culham Science Centre
        effrfss : output real : lower hybrid current drive efficiency (A/W)
        This routine calculates the current drive parameters for a
        lower hybrid system, based on the AEA FUS 172 model.
        AEA FUS 172: Physics Assessment for the European Reactor Study
        """
        rratio = self.lhrad()
        rpenet = rratio * physics_variables.rminor

        # Local density, temperature, toroidal field at this minor radius

        dlocal = 1.0e-19 * self.plasma_profile.neprofile.calculate_profile_y(
            rratio,
            physics_variables.radius_plasma_pedestal_density_norm,
            physics_variables.nd_plasma_electron_on_axis,
            physics_variables.nd_plasma_pedestal_electron,
            physics_variables.nd_plasma_separatrix_electron,
            physics_variables.alphan,
        )
        tlocal = self.plasma_profile.teprofile.calculate_profile_y(
            rratio,
            physics_variables.radius_plasma_pedestal_temp_norm,
            physics_variables.temp_plasma_electron_on_axis_kev,
            physics_variables.temp_plasma_pedestal_kev,
            physics_variables.temp_plasma_separatrix_kev,
            physics_variables.alphat,
            physics_variables.tbeta,
        )
        blocal = (
            physics_variables.b_plasma_toroidal_on_axis
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

        term01 = 6.1e0 / (
            nplacc
            * nplacc
            * (physics_variables.n_charge_plasma_effective_vol_avg + 5.0e0)
        )
        term02 = 1.0e0 + (tlocal / 25.0e0) ** 1.16e0
        term03 = epslh**0.77e0 * np.sqrt(12.25e0 + x * x)
        term04 = 3.5e0 * epslh**0.77e0 + x

        if term03 > term04:
            raise ProcessValueError(
                "Normalised LH efficiency < 0; use a different value of i_hcd_primary",
                term03=term03,
                term04=term04,
            )

        gamlh = term01 * term02 * (1.0e0 - term03 / term04)

        # Current drive efficiency (A/W)

        return gamlh / ((0.1e0 * dlocal) * physics_variables.rmajor)

    def lhrad(self):
        """Routine to calculate Lower Hybrid wave absorption radius
        author: P J Knight, CCFE, Culham Science Centre
        rratio  : output real : minor radius of penetration / rminor
        This routine determines numerically the minor radius at which the
        damping of Lower Hybrid waves occurs, using a Newton-Raphson method.
        AEA FUS 172: Physics Assessment for the European Reactor Study
        """
        #  Correction to refractive index (kept within valid bounds)
        drfind = min(
            0.7e0,
            max(0.1e0, 12.5e0 / physics_variables.temp_plasma_electron_on_axis_kev),
        )

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
            logger.error(
                "LH penetration radius not found after lapno iterations, using 0.8*rminor"
            )
            rat0 = 0.8e0

        return rat0

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
        AEA FUS 172: Physics Assessment for the European Reactor Study
        """
        dlocal = 1.0e-19 * self.plasma_profile.neprofile.calculate_profile_y(
            rratio,
            physics_variables.radius_plasma_pedestal_density_norm,
            physics_variables.nd_plasma_electron_on_axis,
            physics_variables.nd_plasma_pedestal_electron,
            physics_variables.nd_plasma_separatrix_electron,
            physics_variables.alphan,
        )

        #  Local electron temperature

        tlocal = self.plasma_profile.teprofile.calculate_profile_y(
            rratio,
            physics_variables.radius_plasma_pedestal_temp_norm,
            physics_variables.temp_plasma_electron_on_axis_kev,
            physics_variables.temp_plasma_pedestal_kev,
            physics_variables.temp_plasma_separatrix_kev,
            physics_variables.alphat,
            physics_variables.tbeta,
        )

        #  Local toroidal field (evaluated at the inboard region of the flux surface)

        blocal = (
            physics_variables.b_plasma_toroidal_on_axis
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

    def lower_hybrid_fenstermacher(
        self, te: float, rmajor: float, dene20: float
    ) -> float:
        """
        Calculate the lower hybrid frequency using the Fenstermacher formula.
        This function computes the lower hybrid frequency based on the electron
        temperature, major radius, and electron density.

        :param te: Volume averaged electron temperature in keV.
        :type te: float
        :param rmajor: Major radius of the plasma in meters.
        :type rmajor: float
        :param dene20: Volume averaged electron density in units of 10^20 m^-3.
        :type dene20: float

        :return: The calculated absolute current drive efficiency in A/W.
        :rtype: float

        :notes:
            - This forumla was originally in the Oak RidgeSystems Code, attributed to Fenstermacher
              and is used in the AEA FUS 172 report.

        :references:
            - T.C. Hender et al., 'Physics Assessment of the European Reactor Study', AEA FUS 172, 1992.

            - R.L.Reid et al, Oak Ridge Report ORNL/FEDC-87-7, 1988
        """

        return (0.36e0 * (1.0e0 + (te / 25.0e0) ** 1.16e0)) / (rmajor * dene20)

    def lower_hybrid_ehst(
        self, te: float, beta: float, rmajor: float, dene20: float, zeff: float
    ) -> float:
        """
        Calculate the Lower Hybrid current drive efficiency using the Ehst model.

        This function computes the current drive efficiency based on the electron
        temperature, beta, major radius, electron density, and effective charge.

        :param te: Volume averaged electron temperature in keV.
        :type te: float
        :param beta: Plasma beta value (ratio of plasma pressure to magnetic pressure).
        :type beta: float
        :param rmajor: Major radius of the plasma in meters.
        :type rmajor: float
        :param dene20: Volume averaged electron density in units of 10^20 m^-3.
        :type dene20: float
        :param zeff: Plasma effective charge.
        :type zeff: float

        :return: The calculated absolute current drive efficiency in A/W.
        :rtype: float

        :notes:

        :references:
            - Ehst, D.A., and Karney, C.F.F., "Lower Hybrid Current Drive in Tokamaks",
              Nuclear Fusion, 31(10), 1933-1949, 1991.
        """
        return (
            ((te**0.77 * (0.034 + 0.196 * beta)) / (rmajor * dene20))
            * (
                32.0 / (5.0 + zeff)
                + 2.0
                + (12.0 * (6.0 + zeff)) / (5.0 + zeff) / (3.0 + zeff)
                + 3.76 / zeff
            )
            / 12.507
        )


class CurrentDrive:
    def __init__(
        self,
        plasma_profile: PlasmaProfile,
        electron_cyclotron: ElectronCyclotron,
        ion_cyclotron: IonCyclotron,
        lower_hybrid: LowerHybrid,
        neutral_beam: NeutralBeam,
        electron_bernstein: ElectronBernstein,
    ):
        self.outfile = constants.NOUT
        self.plasma_profile = plasma_profile
        self.electron_cyclotron = electron_cyclotron
        self.ion_cyclotron = ion_cyclotron
        self.lower_hybrid = lower_hybrid
        self.neutral_beam = neutral_beam
        self.electron_bernstein = electron_bernstein

    def cudriv(self) -> None:
        """
        Calculate the current drive power requirements.

        This method computes the power requirements of the current drive system
        using a choice of models for the current drive efficiency.

        :param output: Flag indicating whether to write results to the output file.
        :type output: bool

        :raises ProcessValueError: If an invalid current drive switch is encountered.
        """

        current_drive_variables.p_hcd_ecrh_injected_total_mw = 0.0e0
        current_drive_variables.p_hcd_beam_injected_total_mw = 0.0e0
        current_drive_variables.p_hcd_lowhyb_injected_total_mw = 0.0e0
        current_drive_variables.p_hcd_icrh_injected_total_mw = 0.0e0
        current_drive_variables.p_hcd_ebw_injected_total_mw = 0.0e0
        current_drive_variables.c_beam_total = 0.0e0
        current_drive_variables.p_beam_orbit_loss_mw = 0.0e0

        pinjmw1 = 0.0
        p_hcd_primary_ions_mw = 0.0
        p_hcd_primary_electrons_mw = 0.0
        p_hcd_secondary_electrons_mw = 0.0
        p_hcd_secondary_ions_mw = 0.0

        # To stop issues with input file we force
        # zero secondary heating if no injection method
        if current_drive_variables.i_hcd_secondary == 0:
            current_drive_variables.p_hcd_secondary_extra_heat_mw = 0.0

        # i_hcd_calculations |  switch for current drive calculation
        # = 0   |  turned off
        # = 1   |  turned on

        if current_drive_variables.i_hcd_calculations != 0:
            # Put electron density in desired units (10^-20 m-3)
            dene20 = physics_variables.nd_plasma_electrons_vol_avg * 1.0e-20

            # Calculate current drive efficiencies
            # ==============================================================

            # Define a dictionary of lambda functions for current drive efficiency models
            hcd_models = {
                1: lambda: self.lower_hybrid.lower_hybrid_fenstermacher(
                    physics_variables.temp_plasma_electron_vol_avg_kev,
                    physics_variables.rmajor,
                    dene20,
                )
                * current_drive_variables.feffcd,
                2: lambda: self.ion_cyclotron.ion_cyclotron_ipdg89(
                    temp_plasma_electron_density_weighted_kev=physics_variables.temp_plasma_electron_density_weighted_kev,
                    zeff=physics_variables.n_charge_plasma_effective_vol_avg,
                    rmajor=physics_variables.rmajor,
                    dene20=dene20,
                )
                * current_drive_variables.feffcd,
                3: lambda: self.electron_cyclotron.electron_cyclotron_fenstermacher(
                    temp_plasma_electron_density_weighted_kev=physics_variables.temp_plasma_electron_density_weighted_kev,
                    rmajor=physics_variables.rmajor,
                    dene20=dene20,
                    dlamee=physics_variables.dlamee,
                )
                * current_drive_variables.feffcd,
                4: lambda: self.lower_hybrid.lower_hybrid_ehst(
                    te=physics_variables.temp_plasma_electron_vol_avg_kev,
                    beta=physics_variables.beta_total_vol_avg,
                    rmajor=physics_variables.rmajor,
                    dene20=dene20,
                    zeff=physics_variables.n_charge_plasma_effective_vol_avg,
                )
                * current_drive_variables.feffcd,
                5: lambda: (
                    self.neutral_beam.iternb()[0] * current_drive_variables.feffcd
                ),
                6: lambda: self.lower_hybrid.cullhy() * current_drive_variables.feffcd,
                7: lambda: self.electron_cyclotron.culecd()
                * current_drive_variables.feffcd,
                8: lambda: (
                    self.neutral_beam.culnbi()[0] * current_drive_variables.feffcd
                ),
                10: lambda: current_drive_variables.eta_cd_norm_ecrh
                / (dene20 * physics_variables.rmajor),
                12: lambda: self.electron_bernstein.electron_bernstein_freethy(
                    te=physics_variables.temp_plasma_electron_vol_avg_kev,
                    rmajor=physics_variables.rmajor,
                    dene20=dene20,
                    b_plasma_toroidal_on_axis=physics_variables.b_plasma_toroidal_on_axis,
                    n_ecrh_harmonic=current_drive_variables.n_ecrh_harmonic,
                    xi_ebw=current_drive_variables.xi_ebw,
                )
                * current_drive_variables.feffcd,
                13: lambda: self.electron_cyclotron.electron_cyclotron_freethy(
                    te=physics_variables.temp_plasma_electron_vol_avg_kev,
                    zeff=physics_variables.n_charge_plasma_effective_vol_avg,
                    rmajor=physics_variables.rmajor,
                    nd_plasma_electrons_vol_avg=physics_variables.nd_plasma_electrons_vol_avg,
                    b_plasma_toroidal_on_axis=physics_variables.b_plasma_toroidal_on_axis,
                    n_ecrh_harmonic=current_drive_variables.n_ecrh_harmonic,
                    i_ecrh_wave_mode=current_drive_variables.i_ecrh_wave_mode,
                )
                * current_drive_variables.feffcd,
            }

            # Assign outputs for models that return multiple values
            if current_drive_variables.i_hcd_secondary in [5, 8]:
                _, f_p_beam_injected_ions, f_p_beam_shine_through = (
                    self.neutral_beam.iternb()
                    if current_drive_variables.i_hcd_secondary == 5
                    else self.neutral_beam.culnbi()
                )
                current_drive_variables.f_p_beam_injected_ions = f_p_beam_injected_ions
                current_drive_variables.f_p_beam_shine_through = f_p_beam_shine_through

            # Calculate eta_cd_hcd_secondary based on the selected model
            if current_drive_variables.i_hcd_secondary in hcd_models:
                current_drive_variables.eta_cd_hcd_secondary = hcd_models[
                    current_drive_variables.i_hcd_secondary
                ]()
            elif current_drive_variables.i_hcd_secondary != 0:
                raise ProcessValueError(
                    f"Current drive switch is invalid: {current_drive_variables.i_hcd_secondary = }"
                )

            # Calculate eta_cd_hcd_primary based on the selected model
            if current_drive_variables.i_hcd_primary in hcd_models:
                current_drive_variables.eta_cd_hcd_primary = hcd_models[
                    current_drive_variables.i_hcd_primary
                ]()
            else:
                raise ProcessValueError(
                    f"Current drive switch is invalid: {current_drive_variables.i_hcd_primary = }"
                )

            # Calculate the normalised current drive efficieny for the primary heating method
            current_drive_variables.eta_cd_norm_hcd_primary = (
                current_drive_variables.eta_cd_hcd_primary
                * (dene20 * physics_variables.rmajor)
            )
            # Calculate the normalised current drive efficieny for the secondary heating method
            current_drive_variables.eta_cd_norm_hcd_secondary = (
                current_drive_variables.eta_cd_hcd_secondary
                * (dene20 * physics_variables.rmajor)
            )

            # Calculate the driven current for the secondary heating method
            current_drive_variables.c_hcd_secondary_driven = (
                current_drive_variables.eta_cd_hcd_secondary
                * current_drive_variables.p_hcd_secondary_injected_mw
                * 1.0e6
            )
            # Calculate the fraction of the plasma current driven by the secondary heating method
            current_drive_variables.f_c_plasma_hcd_secondary = (
                current_drive_variables.c_hcd_secondary_driven
                / physics_variables.plasma_current
            )

            # Calculate the injected power for the primary heating method
            current_drive_variables.p_hcd_primary_injected_mw = (
                1.0e-6
                * (
                    physics_variables.f_c_plasma_auxiliary
                    - current_drive_variables.f_c_plasma_hcd_secondary
                )
                * physics_variables.plasma_current
                / current_drive_variables.eta_cd_hcd_primary
            )

            # Calculate the driven current for the primary heating method
            current_drive_variables.c_hcd_primary_driven = (
                current_drive_variables.eta_cd_hcd_primary
                * current_drive_variables.p_hcd_primary_injected_mw
                * 1.0e6
            )
            # Calculate the fraction of the plasma current driven by the primary heating method
            current_drive_variables.f_c_plasma_hcd_primary = (
                current_drive_variables.c_hcd_primary_driven
                / physics_variables.plasma_current
            )

            # Calculate the dimensionless current drive efficiency for the primary heating method (ζ)
            current_drive_variables.eta_cd_dimensionless_hcd_primary = self.calculate_dimensionless_current_drive_efficiency(
                nd_plasma_electrons_vol_avg=physics_variables.nd_plasma_electrons_vol_avg,
                rmajor=physics_variables.rmajor,
                temp_plasma_electron_vol_avg_kev=physics_variables.temp_plasma_electron_vol_avg_kev,
                c_hcd_driven=current_drive_variables.c_hcd_primary_driven,
                p_hcd_injected=current_drive_variables.p_hcd_primary_injected_mw
                * 1.0e6,
            )

            if current_drive_variables.p_hcd_secondary_injected_mw > 0.0:
                # Calculate the dimensionless current drive efficiency for the secondary heating method (ζ)
                current_drive_variables.eta_cd_dimensionless_hcd_secondary = self.calculate_dimensionless_current_drive_efficiency(
                    nd_plasma_electrons_vol_avg=physics_variables.nd_plasma_electrons_vol_avg,
                    rmajor=physics_variables.rmajor,
                    temp_plasma_electron_vol_avg_kev=physics_variables.temp_plasma_electron_vol_avg_kev,
                    c_hcd_driven=current_drive_variables.c_hcd_secondary_driven,
                    p_hcd_injected=current_drive_variables.p_hcd_secondary_injected_mw
                    * 1.0e6,
                )

            # ===========================================================

            # Calculate the wall plug power for the secondary heating method
            # ==============================================================

            # Lower hybrid cases
            if current_drive_variables.i_hcd_secondary in [1, 4, 6]:
                # Injected power
                p_hcd_secondary_electrons_mw = (
                    current_drive_variables.p_hcd_secondary_injected_mw
                    + current_drive_variables.p_hcd_secondary_extra_heat_mw
                )

                # Wall plug power
                heat_transport_variables.p_hcd_secondary_electric_mw = (
                    current_drive_variables.p_hcd_secondary_injected_mw
                    + current_drive_variables.p_hcd_secondary_extra_heat_mw
                ) / current_drive_variables.eta_lowhyb_injector_wall_plug

                # Wall plug to injector efficiency
                current_drive_variables.eta_hcd_secondary_injector_wall_plug = (
                    current_drive_variables.eta_lowhyb_injector_wall_plug
                )

                current_drive_variables.p_hcd_lowhyb_injected_total_mw += (
                    current_drive_variables.p_hcd_secondary_injected_mw
                    + current_drive_variables.p_hcd_secondary_extra_heat_mw
                )

            # ==========================================================

            # Ion cyclotron cases
            if current_drive_variables.i_hcd_secondary in [2]:
                # Injected power
                p_hcd_secondary_ions_mw = (
                    current_drive_variables.p_hcd_secondary_injected_mw
                    + current_drive_variables.p_hcd_secondary_extra_heat_mw
                )

                # Wall plug power
                heat_transport_variables.p_hcd_secondary_electric_mw = (
                    current_drive_variables.p_hcd_secondary_injected_mw
                    + current_drive_variables.p_hcd_secondary_extra_heat_mw
                ) / current_drive_variables.eta_icrh_injector_wall_plug

                # Wall plug to injector efficiency
                current_drive_variables.eta_hcd_secondary_injector_wall_plug = (
                    current_drive_variables.eta_icrh_injector_wall_plug
                )

                current_drive_variables.p_hcd_icrh_injected_total_mw += (
                    current_drive_variables.p_hcd_secondary_injected_mw
                    + current_drive_variables.p_hcd_secondary_extra_heat_mw
                )

            # ==========================================================

            # Electron cyclotron cases
            if current_drive_variables.i_hcd_secondary in [3, 7, 10, 13]:
                # Injected power
                p_hcd_secondary_electrons_mw = (
                    current_drive_variables.p_hcd_secondary_injected_mw
                    + current_drive_variables.p_hcd_secondary_extra_heat_mw
                )

                # Wall plug power
                heat_transport_variables.p_hcd_secondary_electric_mw = (
                    current_drive_variables.p_hcd_secondary_injected_mw
                    + current_drive_variables.p_hcd_secondary_extra_heat_mw
                ) / current_drive_variables.eta_ecrh_injector_wall_plug

                # Wall plug to injector efficiency
                current_drive_variables.eta_hcd_secondary_injector_wall_plug = (
                    current_drive_variables.eta_ecrh_injector_wall_plug
                )

                current_drive_variables.p_hcd_ecrh_injected_total_mw += (
                    current_drive_variables.p_hcd_secondary_injected_mw
                    + current_drive_variables.p_hcd_secondary_extra_heat_mw
                )

            # ==========================================================

            # Electron berstein cases
            if current_drive_variables.i_hcd_secondary in [12]:
                # Injected power
                p_hcd_secondary_electrons_mw = (
                    current_drive_variables.p_hcd_secondary_injected_mw
                    + current_drive_variables.p_hcd_secondary_extra_heat_mw
                )

                # Wall plug power
                heat_transport_variables.p_hcd_secondary_electric_mw = (
                    current_drive_variables.p_hcd_secondary_injected_mw
                    + current_drive_variables.p_hcd_secondary_extra_heat_mw
                ) / current_drive_variables.eta_ebw_injector_wall_plug

                # Wall plug to injector efficiency
                current_drive_variables.eta_hcd_secondary_injector_wall_plug = (
                    current_drive_variables.eta_ebw_injector_wall_plug
                )

                current_drive_variables.p_hcd_ebw_injected_total_mw += (
                    current_drive_variables.p_hcd_secondary_injected_mw
                    + current_drive_variables.p_hcd_secondary_extra_heat_mw
                )

            # ==========================================================

            # Neutral beam cases
            elif current_drive_variables.i_hcd_secondary in [5, 8]:
                # Account for first orbit losses
                # (power due to particles that are ionised but not thermalised) [MW]:
                # This includes a second order term in shinethrough*(first orbit loss)

                current_drive_variables.f_p_beam_orbit_loss = min(
                    0.999, current_drive_variables.f_p_beam_orbit_loss
                )  # Should never be needed

                # Shinethrough power (atoms that are not ionised) [MW]:
                current_drive_variables.p_beam_shine_through_mw = (
                    current_drive_variables.p_hcd_secondary_injected_mw
                    + current_drive_variables.p_hcd_secondary_extra_heat_mw
                    * (1.0 - current_drive_variables.f_p_beam_shine_through)
                )

                # First orbit loss
                current_drive_variables.p_beam_orbit_loss_mw = (
                    current_drive_variables.f_p_beam_orbit_loss
                    * (
                        current_drive_variables.p_hcd_secondary_injected_mw
                        + current_drive_variables.p_hcd_secondary_extra_heat_mw
                        - current_drive_variables.p_beam_shine_through_mw
                    )
                )

                # Power deposited
                current_drive_variables.p_beam_plasma_coupled_mw = (
                    current_drive_variables.p_hcd_secondary_injected_mw
                    + current_drive_variables.p_hcd_secondary_extra_heat_mw
                    - current_drive_variables.p_beam_shine_through_mw
                    - current_drive_variables.p_beam_orbit_loss_mw
                )

                p_hcd_secondary_ions_mw = (
                    current_drive_variables.p_beam_plasma_coupled_mw
                    * current_drive_variables.f_p_beam_injected_ions
                )

                p_hcd_secondary_electrons_mw = (
                    current_drive_variables.p_beam_plasma_coupled_mw
                    * (1.0e0 - current_drive_variables.f_p_beam_injected_ions)
                )

                current_drive_variables.pwpnb = (
                    (
                        current_drive_variables.p_hcd_secondary_injected_mw
                        + current_drive_variables.p_hcd_secondary_extra_heat_mw
                    )
                    / current_drive_variables.eta_beam_injector_wall_plug
                )  # neutral beam wall plug power

                heat_transport_variables.p_hcd_secondary_electric_mw = (
                    current_drive_variables.pwpnb
                )

                current_drive_variables.eta_hcd_secondary_injector_wall_plug = (
                    current_drive_variables.eta_beam_injector_wall_plug
                )

                current_drive_variables.c_beam_total = (
                    1.0e-3
                    * (
                        (
                            current_drive_variables.p_hcd_secondary_injected_mw
                            + current_drive_variables.p_hcd_secondary_extra_heat_mw
                        )
                        * 1.0e6
                    )
                    / current_drive_variables.e_beam_kev
                )  # Neutral beam current (A)

                current_drive_variables.p_hcd_beam_injected_total_mw += (
                    current_drive_variables.p_hcd_secondary_injected_mw
                    + current_drive_variables.p_hcd_secondary_extra_heat_mw
                )

            # ==========================================================

            # Lower hybrid cases
            if current_drive_variables.i_hcd_primary in [1, 4, 6]:
                p_hcd_primary_electrons_mw = (
                    current_drive_variables.p_hcd_primary_injected_mw
                    + current_drive_variables.p_hcd_primary_extra_heat_mw
                )

                current_drive_variables.p_hcd_lowhyb_injected_total_mw += (
                    current_drive_variables.p_hcd_primary_injected_mw
                    + current_drive_variables.p_hcd_primary_extra_heat_mw
                )

                # Wall plug power
                heat_transport_variables.p_hcd_primary_electric_mw = (
                    current_drive_variables.p_hcd_primary_injected_mw
                    + current_drive_variables.p_hcd_primary_extra_heat_mw
                ) / current_drive_variables.eta_lowhyb_injector_wall_plug

                # Wall plug to injector efficiency
                current_drive_variables.eta_hcd_primary_injector_wall_plug = (
                    current_drive_variables.eta_lowhyb_injector_wall_plug
                )

                # Wall plug power
                current_drive_variables.p_hcd_lowhyb_electric_mw = (
                    current_drive_variables.p_hcd_lowhyb_injected_total_mw
                    / current_drive_variables.eta_lowhyb_injector_wall_plug
                )

            # ===========================================================

            # Ion cyclotron cases
            if current_drive_variables.i_hcd_primary in [2]:
                p_hcd_primary_ions_mw = (
                    current_drive_variables.p_hcd_primary_injected_mw
                    + current_drive_variables.p_hcd_primary_extra_heat_mw
                )

                # Wall plug power
                heat_transport_variables.p_hcd_primary_electric_mw = (
                    current_drive_variables.p_hcd_primary_injected_mw
                    + current_drive_variables.p_hcd_primary_extra_heat_mw
                ) / current_drive_variables.eta_icrh_injector_wall_plug

                # Wall plug to injector efficiency
                current_drive_variables.eta_hcd_primary_injector_wall_plug = (
                    current_drive_variables.eta_icrh_injector_wall_plug
                )

                current_drive_variables.p_hcd_icrh_injected_total_mw += (
                    current_drive_variables.p_hcd_primary_injected_mw
                    + current_drive_variables.p_hcd_primary_extra_heat_mw
                )

                # Wall plug power
                current_drive_variables.p_hcd_icrh_electric_mw = (
                    current_drive_variables.p_hcd_icrh_injected_total_mw
                    / current_drive_variables.eta_icrh_injector_wall_plug
                )

            # ===========================================================

            # Electron cyclotron cases

            if current_drive_variables.i_hcd_primary in [3, 7, 10, 13]:
                p_hcd_primary_electrons_mw = (
                    current_drive_variables.p_hcd_primary_injected_mw
                    + current_drive_variables.p_hcd_primary_extra_heat_mw
                )

                # Wall plug to injector efficiency
                heat_transport_variables.p_hcd_primary_electric_mw = (
                    current_drive_variables.p_hcd_primary_injected_mw
                    + current_drive_variables.p_hcd_primary_extra_heat_mw
                ) / current_drive_variables.eta_ecrh_injector_wall_plug

                current_drive_variables.eta_hcd_primary_injector_wall_plug = (
                    current_drive_variables.eta_ecrh_injector_wall_plug
                )

                current_drive_variables.p_hcd_ecrh_injected_total_mw += (
                    current_drive_variables.p_hcd_primary_injected_mw
                    + current_drive_variables.p_hcd_primary_extra_heat_mw
                )

                # Wall plug power
                current_drive_variables.p_hcd_ecrh_electric_mw = (
                    current_drive_variables.p_hcd_ecrh_injected_total_mw
                    / current_drive_variables.eta_ecrh_injector_wall_plug
                )

            # ===========================================================

            # Electron bernstein cases

            if current_drive_variables.i_hcd_primary in [12]:
                p_hcd_primary_electrons_mw = (
                    current_drive_variables.p_ebw_injected_mw
                    + current_drive_variables.p_hcd_primary_extra_heat_mw
                )

                # Wall plug to injector efficiency
                heat_transport_variables.p_hcd_primary_electric_mw = (
                    current_drive_variables.p_hcd_primary_injected_mw
                    + current_drive_variables.p_hcd_primary_extra_heat_mw
                ) / current_drive_variables.eta_ebw_injector_wall_plug

                current_drive_variables.eta_hcd_primary_injector_wall_plug = (
                    current_drive_variables.eta_ebw_injector_wall_plug
                )

                current_drive_variables.p_hcd_ebw_injected_total_mw += (
                    current_drive_variables.p_hcd_primary_injected_mw
                    + current_drive_variables.p_hcd_primary_extra_heat_mw
                )

                # Wall plug power
                current_drive_variables.p_hcd_ebw_electric_mw = (
                    current_drive_variables.p_ebw_injected_mw
                    / current_drive_variables.eta_ebw_injector_wall_plug
                )

            # ===========================================================

            elif current_drive_variables.i_hcd_primary in [5, 8]:
                # Account for first orbit losses
                # (power due to particles that are ionised but not thermalised) [MW]:
                # This includes a second order term in shinethrough*(first orbit loss)
                current_drive_variables.f_p_beam_orbit_loss = min(
                    0.999, current_drive_variables.f_p_beam_orbit_loss
                )  # Should never be needed

                # Shinethrough power (atoms that are not ionised) [MW]:
                current_drive_variables.p_beam_shine_through_mw = (
                    current_drive_variables.p_hcd_primary_injected_mw
                    + current_drive_variables.p_hcd_primary_extra_heat_mw
                    * (1.0 - current_drive_variables.f_p_beam_shine_through)
                )

                # First orbit loss
                current_drive_variables.p_beam_orbit_loss_mw = (
                    current_drive_variables.f_p_beam_orbit_loss
                    * (
                        current_drive_variables.p_hcd_primary_injected_mw
                        + current_drive_variables.p_hcd_primary_extra_heat_mw
                        - current_drive_variables.p_beam_shine_through_mw
                    )
                )

                # Power deposited
                current_drive_variables.p_beam_plasma_coupled_mw = (
                    current_drive_variables.p_hcd_primary_injected_mw
                    + current_drive_variables.p_hcd_primary_extra_heat_mw
                    - current_drive_variables.p_beam_shine_through_mw
                    - current_drive_variables.p_beam_orbit_loss_mw
                )

                p_hcd_primary_ions_mw = (
                    pinjmw1 * current_drive_variables.f_p_beam_injected_ions
                )
                p_hcd_primary_electrons_mw = pinjmw1 * (
                    1.0e0 - current_drive_variables.f_p_beam_injected_ions
                )

                current_drive_variables.pwpnb = (
                    current_drive_variables.p_hcd_primary_injected_mw
                    + current_drive_variables.p_hcd_primary_extra_heat_mw
                    / current_drive_variables.eta_beam_injector_wall_plug
                )

                # Neutral beam wall plug power
                heat_transport_variables.p_hcd_primary_electric_mw = (
                    current_drive_variables.pwpnb
                )
                current_drive_variables.eta_hcd_primary_injector_wall_plug = (
                    current_drive_variables.eta_beam_injector_wall_plug
                )

                current_drive_variables.c_beam_total = (
                    1.0e-3
                    * (
                        (
                            current_drive_variables.p_hcd_primary_injected_mw
                            + current_drive_variables.p_hcd_primary_extra_heat_mw
                        )
                        * 1.0e6
                    )
                    / current_drive_variables.e_beam_kev
                )  # Neutral beam current (A)

                current_drive_variables.p_hcd_beam_injected_total_mw += (
                    current_drive_variables.p_hcd_primary_injected_mw
                    + current_drive_variables.p_hcd_primary_extra_heat_mw
                )

            # ===========================================================

            # Total injected power that contributed to heating
            current_drive_variables.p_hcd_injected_total_mw = (
                current_drive_variables.p_hcd_primary_injected_mw
                + current_drive_variables.p_hcd_primary_extra_heat_mw
                + current_drive_variables.p_hcd_secondary_injected_mw
                + current_drive_variables.p_hcd_secondary_extra_heat_mw
            )

            # Total injected power that contributed to current drive
            current_drive_variables.p_hcd_injected_current_total_mw = (
                current_drive_variables.p_hcd_primary_injected_mw
                + current_drive_variables.p_hcd_secondary_injected_mw
            )

            pinjmw1 = p_hcd_primary_electrons_mw + p_hcd_primary_ions_mw

            # Total injected power given to electrons
            current_drive_variables.p_hcd_injected_electrons_mw = (
                p_hcd_primary_electrons_mw + p_hcd_secondary_electrons_mw
            )

            # Total injected power given to ions
            current_drive_variables.p_hcd_injected_ions_mw = (
                p_hcd_primary_ions_mw + p_hcd_secondary_ions_mw
            )

            # Total wall plug power for all heating systems
            heat_transport_variables.p_hcd_electric_total_mw = (
                heat_transport_variables.p_hcd_primary_electric_mw
                + heat_transport_variables.p_hcd_secondary_electric_mw
            )

            # Reset injected power to zero for ignited plasma (fudge)
            if physics_variables.i_plasma_ignited == 1:
                heat_transport_variables.p_hcd_electric_total_mw = 0.0e0

            # Ratio of fusion to input (injection+ohmic) power
            current_drive_variables.big_q_plasma = (
                physics_variables.p_fusion_total_mw
                / (
                    current_drive_variables.p_hcd_injected_total_mw
                    + current_drive_variables.p_beam_orbit_loss_mw
                    + physics_variables.p_plasma_ohmic_mw
                )
            )

    def calculate_dimensionless_current_drive_efficiency(
        self,
        nd_plasma_electrons_vol_avg: float,
        rmajor: float,
        temp_plasma_electron_vol_avg_kev: float,
        c_hcd_driven: float,
        p_hcd_injected: float,
    ) -> float:
        """
                Calculate the dimensionless current drive efficiency, ζ.

                This function computes the dimensionless current drive efficiency
                based on the average electron density, major radius, and electron temperature.

                :param nd_plasma_electrons_vol_avg: Volume averaged electron density in m^-3.
                :type nd_plasma_electrons_vol_avg: float
                :param rmajor: Major radius of the plasma in meters.
                :type rmajor: float
                :param temp_plasma_electron_vol_avg_kev: Volume averaged electron temperature in keV.
                :type temp_plasma_electron_vol_avg_kev: float
                :param c_hcd_driven: Current driven by the heating and current drive system.
                :type c_hcd_driven: float
                :param p_hcd_injected: Power injected by the heating and current drive system.
                :type p_hcd_injected: float
                :return: The calculated dimensionless current drive efficiency.
                :rtype: float

                :references:
                    - E. Poli et al., “Electron-cyclotron-current-drive efficiency in DEMO plasmas,”
                    Nuclear Fusion, vol. 53, no. 1, pp. 013011-013011, Dec. 2012,
                    doi: https://doi.org/10.1088/0029-5515/53/1/013011.
        ‌
                    - T. C. Luce et al., “Generation of Localized Noninductive Current by Electron Cyclotron Waves on the DIII-D Tokamak,”
                    Physical Review Letters, vol. 83, no. 22, pp. 4550-4553, Nov. 1999,
                    doi: https://doi.org/10.1103/physrevlett.83.4550.
        """

        return (
            (constants.ELECTRON_CHARGE**3 / constants.EPSILON0**2)
            * (
                (nd_plasma_electrons_vol_avg * rmajor)
                / (temp_plasma_electron_vol_avg_kev * constants.KILOELECTRON_VOLT)
            )
            * (c_hcd_driven / p_hcd_injected)
        )

    def output_current_drive(self):
        """
        Output the current drive information to the output file.
        This method writes the current drive information to the output file.
        """

        po.oheadr(self.outfile, "Heating & Current Drive System")

        if physics_variables.i_plasma_ignited == 1:
            po.ocmmnt(
                self.outfile,
                "Ignited plasma; injected power only used for start-up phase",
            )

        if abs(physics_variables.f_c_plasma_inductive) > 1.0e-8:
            po.ocmmnt(
                self.outfile,
                "Current is driven by both inductive and non-inductive means.",
            )
            po.oblnkl(self.outfile)

        po.ovarre(
            self.outfile,
            "Fusion gain factor Q",
            "(big_q_plasma)",
            current_drive_variables.big_q_plasma,
            "OP ",
        )
        po.oblnkl(self.outfile)

        if current_drive_variables.i_hcd_calculations == 0:
            po.ocmmnt(self.outfile, "No current drive used")
            po.oblnkl(self.outfile)
            return

        po.ovarin(
            self.outfile,
            "Primary current drive efficiency model",
            "(i_hcd_primary)",
            current_drive_variables.i_hcd_primary,
        )

        if current_drive_variables.i_hcd_primary in [1, 4, 6]:
            po.ocmmnt(self.outfile, "Lower Hybrid Current Drive")
        elif current_drive_variables.i_hcd_primary == 2:
            po.ocmmnt(self.outfile, "Ion Cyclotron Current Drive")
        elif current_drive_variables.i_hcd_primary in [3, 7]:
            po.ocmmnt(self.outfile, "Electron Cyclotron Current Drive")
        elif current_drive_variables.i_hcd_primary in [5, 8]:
            po.ocmmnt(self.outfile, "Neutral Beam Current Drive")
        elif current_drive_variables.i_hcd_primary == 10:
            po.ocmmnt(
                self.outfile,
                "Electron Cyclotron Current Drive (input normalised efficiency)",
            )
        elif current_drive_variables.i_hcd_primary == 12:
            po.ocmmnt(self.outfile, "Electron Bernstein Wave Current Drive")
        elif current_drive_variables.i_hcd_primary == 13:
            po.ocmmnt(
                self.outfile,
                "Electron Cyclotron Current Drive (with Zeff & Te dependance)",
            )

        po.oblnkl(self.outfile)

        po.ovarre(
            self.outfile,
            "Absolute current drive efficiency of primary system [A/W]",
            "(eta_cd_hcd_primary)",
            current_drive_variables.eta_cd_hcd_primary,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Normalised current drive efficiency of primary system [10^20 A / Wm^2]",
            "(eta_cd_norm_hcd_primary)",
            current_drive_variables.eta_cd_norm_hcd_primary,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Dimensionless current drive efficiency of primary system, ζ",
            "(eta_cd_dimensionless_hcd_primary)",
            current_drive_variables.eta_cd_dimensionless_hcd_primary,
            "OP ",
        )
        if current_drive_variables.i_hcd_primary == 10:
            po.ovarre(
                self.outfile,
                "ECRH plasma heating efficiency",
                "(eta_cd_norm_ecrh)",
                current_drive_variables.eta_cd_norm_ecrh,
            )
        po.ovarre(
            self.outfile,
            "Power injected into plasma by primary system for current drive (MW)",
            "(p_hcd_primary_injected_mw)",
            current_drive_variables.p_hcd_primary_injected_mw,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Extra power injected into plasma by primary system  (MW)",
            "(p_hcd_primary_extra_heat_mw)",
            current_drive_variables.p_hcd_primary_extra_heat_mw,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Current driven in plasma by primary system (A)",
            "(c_hcd_primary_driven)",
            current_drive_variables.c_hcd_primary_driven,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Fraction of plasma current driven by primary system",
            "(f_c_plasma_hcd_primary)",
            current_drive_variables.f_c_plasma_hcd_primary,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Wall plug to injector efficiency of primary system",
            "(eta_hcd_primary_injector_wall_plug)",
            current_drive_variables.eta_hcd_primary_injector_wall_plug,
            "IP ",
        )
        po.ovarre(
            self.outfile,
            "Wall plug electric power of primary system",
            "(p_hcd_primary_electric_mw)",
            heat_transport_variables.p_hcd_primary_electric_mw,
            "OP ",
        )

        if current_drive_variables.i_hcd_primary in [12, 13]:
            po.oblnkl(self.outfile)
            po.ovarre(
                self.outfile,
                "ECRH / EBW harmonic number",
                "(n_ecrh_harmonic)",
                current_drive_variables.n_ecrh_harmonic,
            )
            po.ovarre(
                self.outfile,
                "EBW coupling efficiency",
                "(xi_ebw)",
                current_drive_variables.xi_ebw,
            )
        if current_drive_variables.i_hcd_primary == 13:
            po.ovarin(
                self.outfile,
                "Electron cyclotron cutoff wave mode switch",
                "(i_ecrh_wave_mode)",
                current_drive_variables.i_ecrh_wave_mode,
            )

        po.oblnkl(self.outfile)

        if current_drive_variables.i_hcd_primary in [5, 8]:
            po.oblnkl(self.outfile)
            po.ocmmnt(self.outfile, "Neutral beam power balance :")
            po.ocmmnt(self.outfile, "----------------------------")

            po.ovarre(
                self.outfile,
                "Neutral beam energy (keV)",
                "(e_beam_kev)",
                current_drive_variables.e_beam_kev,
            )
            if (current_drive_variables.i_hcd_primary == 5) or (
                current_drive_variables.i_hcd_primary == 8
            ):
                po.ovarre(
                    self.outfile,
                    "Neutral beam current (A)",
                    "(c_beam_total)",
                    current_drive_variables.c_beam_total,
                    "OP ",
                )

            po.ovarre(
                self.outfile,
                "Neutral beam wall plug efficiency",
                "(eta_beam_injector_wall_plug)",
                current_drive_variables.eta_beam_injector_wall_plug,
            )
            po.ovarre(
                self.outfile,
                "Beam decay lengths to centre",
                "(n_beam_decay_lengths_core)",
                current_drive_variables.n_beam_decay_lengths_core,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Beam shine-through fraction",
                "(f_p_beam_shine_through)",
                current_drive_variables.f_p_beam_shine_through,
                "OP ",
            )

            if (current_drive_variables.i_hcd_primary == 5) or (
                current_drive_variables.i_hcd_primary == 8
            ):
                po.ovarrf(
                    self.outfile,
                    "Beam first orbit loss power (MW)",
                    "(p_beam_orbit_loss_mw)",
                    current_drive_variables.p_beam_orbit_loss_mw,
                    "OP ",
                )
                po.ovarrf(
                    self.outfile,
                    "Beam shine-through power [MW]",
                    "(p_beam_shine_through_mw)",
                    current_drive_variables.p_beam_shine_through_mw,
                    "OP ",
                )
                po.ovarrf(
                    self.outfile,
                    "Maximum allowable beam power (MW)",
                    "(p_hcd_injected_max)",
                    current_drive_variables.p_hcd_injected_max,
                )
                po.oblnkl(self.outfile)
                po.ovarrf(
                    self.outfile,
                    "Beam power entering vacuum vessel (MW)",
                    "(p_beam_injected_mw)",
                    current_drive_variables.p_beam_injected_mw,
                    "OP ",
                )
                po.ovarre(
                    self.outfile,
                    "Fraction of beam energy to ions",
                    "(f_p_beam_injected_ions)",
                    current_drive_variables.f_p_beam_injected_ions,
                    "OP ",
                )
                po.ovarre(
                    self.outfile,
                    "Beam duct shielding thickness (m)",
                    "(dx_beam_shield)",
                    current_drive_variables.dx_beam_shield,
                )
                po.ovarre(
                    self.outfile,
                    "Beam tangency radius / Plasma major radius",
                    "(f_radius_beam_tangency_rmajor)",
                    current_drive_variables.f_radius_beam_tangency_rmajor,
                )
                po.ovarre(
                    self.outfile,
                    "Beam centreline tangency radius (m)",
                    "(radius_beam_tangency)",
                    current_drive_variables.radius_beam_tangency,
                    "OP ",
                )
                po.ovarre(
                    self.outfile,
                    "Maximum possible tangency radius (m)",
                    "(radius_beam_tangency_max)",
                    current_drive_variables.radius_beam_tangency_max,
                    "OP ",
                )

        po.ocmmnt(self.outfile, "----------------------------")
        po.oblnkl(self.outfile)
        po.ovarin(
            self.outfile,
            "Secondary current drive efficiency model",
            "(i_hcd_secondary)",
            current_drive_variables.i_hcd_secondary,
        )

        if current_drive_variables.i_hcd_secondary in [1, 4, 6]:
            po.ocmmnt(self.outfile, "Lower Hybrid Current Drive")
        elif current_drive_variables.i_hcd_secondary == 2:
            po.ocmmnt(self.outfile, "Ion Cyclotron Current Drive")
        elif current_drive_variables.i_hcd_secondary in [3, 7]:
            po.ocmmnt(self.outfile, "Electron Cyclotron Current Drive")
        elif current_drive_variables.i_hcd_secondary in [5, 8]:
            po.ocmmnt(self.outfile, "Neutral Beam Current Drive")
        elif current_drive_variables.i_hcd_secondary == 10:
            po.ocmmnt(
                self.outfile,
                "Electron Cyclotron Current Drive (input normalised efficiency)",
            )
        elif current_drive_variables.i_hcd_secondary == 12:
            po.ocmmnt(self.outfile, "Electron Bernstein Wave Current Drive")
        elif current_drive_variables.i_hcd_secondary == 13:
            po.ocmmnt(
                self.outfile,
                "Electron Cyclotron Current Drive (with Zeff & Te dependance)",
            )
        po.oblnkl(self.outfile)

        po.ovarre(
            self.outfile,
            "Absolute current drive efficiency of secondary system [A/W]",
            "(eta_cd_hcd_secondary)",
            current_drive_variables.eta_cd_hcd_secondary,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Normalised current drive efficiency of secondary system [10^20 A / Wm^2]",
            "(eta_cd_norm_hcd_secondary)",
            current_drive_variables.eta_cd_norm_hcd_secondary,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Dimensionless current drive efficiency of secondary system, ζ",
            "(eta_cd_dimensionless_hcd_secondary)",
            current_drive_variables.eta_cd_dimensionless_hcd_secondary,
            "OP ",
        )
        if current_drive_variables.i_hcd_secondary == 10:
            po.ovarre(
                self.outfile,
                "ECRH plasma heating efficiency",
                "(eta_cd_norm_ecrh)",
                current_drive_variables.eta_cd_norm_ecrh,
            )

        po.ovarre(
            self.outfile,
            "Power injected into plasma by secondary system (MW)",
            "(p_hcd_secondary_injected_mw)",
            current_drive_variables.p_hcd_secondary_injected_mw,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Extra power injected into plasma by secondary system  (MW)",
            "(p_hcd_secondary_extra_heat_mw)",
            current_drive_variables.p_hcd_secondary_extra_heat_mw,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Current driven in plasma by secondary system (A)",
            "(c_hcd_secondary_driven)",
            current_drive_variables.c_hcd_secondary_driven,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Fraction of plasma current driven by secondary system",
            "(f_c_plasma_hcd_secondary)",
            current_drive_variables.f_c_plasma_hcd_secondary,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Wall plug to injector efficiency of secondary system",
            "(eta_hcd_secondary_injector_wall_plug)",
            current_drive_variables.eta_hcd_secondary_injector_wall_plug,
            "IP ",
        )
        po.ovarre(
            self.outfile,
            "Wall plug electric power of secondary system",
            "(p_hcd_secondary_electric_mw)",
            heat_transport_variables.p_hcd_secondary_electric_mw,
            "OP ",
        )

        po.oblnkl(self.outfile)

        if current_drive_variables.i_hcd_secondary in [5, 8]:
            po.oblnkl(self.outfile)
            po.ocmmnt(self.outfile, "Neutral beam power balance :")
            po.ocmmnt(self.outfile, "----------------------------")

            po.ovarre(
                self.outfile,
                "Neutral beam energy (keV)",
                "(e_beam_kev)",
                current_drive_variables.e_beam_kev,
            )
            if (current_drive_variables.i_hcd_primary == 5) or (
                current_drive_variables.i_hcd_primary == 8
            ):
                po.ovarre(
                    self.outfile,
                    "Neutral beam current (A)",
                    "(c_beam_total)",
                    current_drive_variables.c_beam_total,
                    "OP ",
                )

            po.ovarre(
                self.outfile,
                "Neutral beam wall plug efficiency",
                "(eta_beam_injector_wall_plug)",
                current_drive_variables.eta_beam_injector_wall_plug,
            )
            po.ovarre(
                self.outfile,
                "Beam decay lengths to centre",
                "(n_beam_decay_lengths_core)",
                current_drive_variables.n_beam_decay_lengths_core,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Beam shine-through fraction",
                "(f_p_beam_shine_through)",
                current_drive_variables.f_p_beam_shine_through,
                "OP ",
            )

            if (current_drive_variables.i_hcd_primary == 5) or (
                current_drive_variables.i_hcd_primary == 8
            ):
                po.ovarrf(
                    self.outfile,
                    "Beam first orbit loss power (MW)",
                    "(p_beam_orbit_loss_mw)",
                    current_drive_variables.p_beam_orbit_loss_mw,
                    "OP ",
                )
                po.ovarrf(
                    self.outfile,
                    "Beam shine-through power [MW]",
                    "(p_beam_shine_through_mw)",
                    current_drive_variables.p_beam_shine_through_mw,
                    "OP ",
                )
                po.ovarrf(
                    self.outfile,
                    "Maximum allowable beam power (MW)",
                    "(p_hcd_injected_max)",
                    current_drive_variables.p_hcd_injected_max,
                )
                po.oblnkl(self.outfile)
                po.ovarrf(
                    self.outfile,
                    "Beam power entering vacuum vessel (MW)",
                    "(p_beam_injected_mw)",
                    current_drive_variables.p_beam_injected_mw,
                    "OP ",
                )
                po.ovarre(
                    self.outfile,
                    "Fraction of beam energy to ions",
                    "(f_p_beam_injected_ions)",
                    current_drive_variables.f_p_beam_injected_ions,
                    "OP ",
                )
                po.ovarre(
                    self.outfile,
                    "Beam duct shielding thickness (m)",
                    "(dx_beam_shield)",
                    current_drive_variables.dx_beam_shield,
                )
                po.ovarre(
                    self.outfile,
                    "Beam tangency radius / Plasma major radius",
                    "(f_radius_beam_tangency_rmajor)",
                    current_drive_variables.f_radius_beam_tangency_rmajor,
                )
                po.ovarre(
                    self.outfile,
                    "Beam centreline tangency radius (m)",
                    "(radius_beam_tangency)",
                    current_drive_variables.radius_beam_tangency,
                    "OP ",
                )
                po.ovarre(
                    self.outfile,
                    "Maximum possible tangency radius (m)",
                    "(radius_beam_tangency_max)",
                    current_drive_variables.radius_beam_tangency_max,
                    "OP ",
                )

        po.ocmmnt(self.outfile, "----------------------------")

        po.osubhd(self.outfile, "Totals :")

        po.ovarre(
            self.outfile,
            "Total injected heating power that drove plasma current (MW)",
            "(p_hcd_injected_current_total_mw)",
            current_drive_variables.p_hcd_injected_current_total_mw,
        )
        po.ovarre(
            self.outfile,
            "Total injected heating power across all systems (MW)",
            "(p_hcd_injected_total_mw)",
            current_drive_variables.p_hcd_injected_total_mw,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Total injected heating power given to the electrons (MW)",
            "(p_hcd_injected_electrons_mw)",
            current_drive_variables.p_hcd_injected_electrons_mw,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Total injected heating power given to the ions (MW)",
            "(p_hcd_injected_ions_mw)",
            current_drive_variables.p_hcd_injected_ions_mw,
            "OP ",
        )

        po.ovarre(
            self.outfile,
            "Upper limit on total plasma injected power (MW)",
            "(p_hcd_injected_max)",
            current_drive_variables.p_hcd_injected_max,
            "OP ",
        )

        po.osubhd(self.outfile, "Contributions:")

        po.ovarre(
            self.outfile,
            "Injected power into plasma from lower hybrid systems (MW)",
            "(p_hcd_lowhyb_injected_total_mw)",
            current_drive_variables.p_hcd_lowhyb_injected_total_mw,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Injected power into plasma from ion cyclotron systems (MW)",
            "(p_hcd_icrh_injected_total_mw)",
            current_drive_variables.p_hcd_icrh_injected_total_mw,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Injected power into plasma from electron cyclotron systems (MW)",
            "(p_hcd_ecrh_injected_total_mw)",
            current_drive_variables.p_hcd_ecrh_injected_total_mw,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Injected power into plasma from neutral beam systems (MW)",
            "(p_hcd_beam_injected_total_mw)",
            current_drive_variables.p_hcd_beam_injected_total_mw,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Injected power into plasma from lower hybrid systems (MW)",
            "(p_hcd_ebw_injected_total_mw)",
            current_drive_variables.p_hcd_ebw_injected_total_mw,
            "OP ",
        )

        po.osubhd(self.outfile, "Fractions of current drive :")
        po.ovarrf(
            self.outfile,
            "Bootstrap fraction",
            "(f_c_plasma_bootstrap)",
            current_drive_variables.f_c_plasma_bootstrap,
            "OP ",
        )
        po.ovarrf(
            self.outfile,
            "Diamagnetic fraction",
            "(f_c_plasma_diamagnetic)",
            current_drive_variables.f_c_plasma_diamagnetic,
            "OP ",
        )
        po.ovarrf(
            self.outfile,
            "Pfirsch-Schlueter fraction",
            "(f_c_plasma_pfirsch_schluter)",
            current_drive_variables.f_c_plasma_pfirsch_schluter,
            "OP ",
        )
        po.ovarrf(
            self.outfile,
            "Auxiliary current drive fraction",
            "(f_c_plasma_auxiliary)",
            physics_variables.f_c_plasma_auxiliary,
            "OP ",
        )
        po.ovarrf(
            self.outfile,
            "Inductive fraction",
            "(f_c_plasma_inductive)",
            physics_variables.f_c_plasma_inductive,
            "OP ",
        )

        # MDK Add physics_variables.f_c_plasma_non_inductive as it can be an iteration variable
        po.ovarrf(
            self.outfile,
            "Fraction of the plasma current produced by non-inductive means",
            "(f_c_plasma_non_inductive)",
            physics_variables.f_c_plasma_non_inductive,
        )

        if (
            abs(
                current_drive_variables.f_c_plasma_bootstrap
                - current_drive_variables.f_c_plasma_bootstrap_max
            )
            < 1.0e-8
        ):
            po.ocmmnt(self.outfile, "Warning : bootstrap current fraction is at")
            po.ocmmnt(self.outfile, "          its prescribed maximum.")

        po.oblnkl(self.outfile)
