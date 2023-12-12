import math

from process.fortran import constants
from process.fortran import build_variables as bv
from process.fortran import divertor_variables as dv
from process.fortran import physics_variables as pv
from process.fortran import tfcoil_variables as tfv
from process.fortran import process_output as po
from process.fortran import error_handling as eh


class Divertor:
    """Module containing divertor routines
    author: P J Knight, CCFE, Culham Science Centre

    This module contains routines relevant for calculating the
    divertor parameters for a fusion power plant.
    AEA FUS 251: A User's Guide to the PROCESS Systems Code
    """

    def __init__(self) -> None:
        self.outfile = constants.nout  # output file unit

    def run(self, output: bool) -> None:
        """Routine to call the divertor model
        author: J Galambos, ORNL
        author: P J Knight, CCFE, Culham Science Centre

        This subroutine calls the divertor routine. This routine scales
        dimensions, powers and field levels which are used as input to
        the Harrison divertor model.
        AEA FUS 251: A User's Guide to the PROCESS Systems Code

        :param output: indicate whether output should be written to the output file, or not
        :type output: boolean
        """

        if pv.itart == 1 and (dv.i_hldiv == 0 or dv.i_hldiv == 1):
            self.divtart(
                pv.rmajor,
                pv.rminor,
                pv.triang,
                bv.scrapli,
                bv.vgap,
                pv.pdivt,
                output=output,
            )
            return
        elif dv.i_hldiv == 2:
            self.divwade(
                pv.rmajor,
                pv.rminor,
                pv.aspect,
                pv.bt,
                pv.bp,
                pv.pdivt,
                dv.flux_exp,
                pv.nesep,
                dv.beta_div,
                pv.rad_fraction_sol,
                pv.ftar,
                output=output,
            )
            return

        #  Scale geometric quantities

        #  Perpendicular diffusivity in the plasma scrapeoff (m2/s)
        #  Assume ions transport 33% of electron power

        xperp = dv.xpertin * 1.33e0

        #  Reference null to strike distances
        #  Only set up for outer divertor for double-null

        # plsep = min(plsepo,pi*rminor) #  Obscure reason to set a limit...
        plsep = bv.plsepo

        #  Scale plasma quantities

        delne = dv.prn1 * pv.dene * 1.0e-20  # scrapeoff density by main plasma
        pwr = pv.pdivt  # power flow to divertor (MW)
        aionso = pv.afuel  # scrape-off layer ion mass

        if dv.divdum == 0:  # Divertor Zeff: scaled
            zeffso = 1.0e0 + 0.8e0 * (pv.zeff - 1.0e0)
        else:  # use input value
            zeffso = dv.zeffdiv

        #  Strike point field values

        bpstk = pv.bp * 0.45e0
        btstk = pv.bt * pv.rmajor / bv.rspo
        rbpbt = bpstk / btstk

        #  Parallel diffusivity in the plasma scrapeoff (m2/s)

        xpara = dv.xparain / zeffso

        #  Null radius

        rnull = pv.rmajor - pv.rminor * pv.triang

        #  Divertor area and radius ratio

        rsrd = (rnull + pv.rmajor + pv.rminor) / (rnull + bv.rspo)
        diva = constants.pi * (rnull + bv.rspo) * plsep
        adas = diva / pv.sarea

        #  Main plasma separatrix area to divertor (and power fraction)
        # +PJK Is the 2 related to 2 divertors (i.e. double-null assumed)?
        frgd = (pv.sareao) / (2.0e0 * pv.sarea)
        # -PJK
        #  Power flow to divertor

        pdiv = pwr * dv.ksic / 2.0e0
        qdiv = pdiv / (pv.sarea * frgd)

        #  Connection length scalings
        #  (2.5 factor comes from normalization to ITER 1990)

        tconl = (
            2.5e0 * pv.rmajor * pv.q * (1.0e0 + 1.0e0 / (pv.q * pv.aspect) ** 2) ** 0.5
        )
        dtheta = plsep / pv.rminor
        dconl = (
            2.5e0
            * bv.rspo
            * pv.q
            * dtheta
            * (1.0e0 + 1.0e0 / (pv.q * pv.aspect) ** 2) ** 0.5
        )
        rconl = dconl / tconl

        #  Minimum strike angle

        minstang = 0.5e0

        #  Call divertor routine

        (
            delta,
            delw,
            dv.dendiv,
            dv.densin,
            gamdt,
            dv.lamp,
            omlarg,
            dv.ppdiv,
            dv.ppdivr,
            dv.ptpdiv,
            dv.tdiv,
            dv.tsep,
        ) = self.divert(
            adas,
            aionso,
            dv.anginc,
            delne,
            dv.c1div,
            dv.c2div,
            dv.c3div,
            dv.c4div,
            dv.c5div,
            dv.delld,
            dv.fdfs,
            dv.fififi,
            frgd,
            dv.frrp,
            minstang,
            dv.omegan,
            qdiv,
            pdiv,
            rbpbt,
            rconl,
            pv.rmajor,
            rsrd,
            tconl,
            xpara,
            xperp,
        )

        #  Heat load

        dv.hldiv = dv.ppdivr

        #  Ratio of collision length to connection length

        dv.rlclolcn = 1.44e-3 * dv.tsep**2 / (delne * 15.0e0 * tconl)

        # output deleted as per discussion in #242.
        # tnunn: 19/01/2022

    def divert(
        self,
        adas,
        aion,
        anginc,
        delne,
        c1div,
        c2div,
        c3div,
        c4div,
        c5div,
        delld,
        fdfs,
        fififi,
        frgd,
        frrp,
        minstang,
        omegan,
        qdiv,
        pdiv,
        rbpbt,
        rconl,
        rmaj,
        rsrd,
        tconl,
        xpara,
        xperp,
    ):
        """Harrison-Kukushkin analytic ITER divertor model
        author: J Galambos, ORNL
        author: P J Knight, CCFE, Culham Science Centre

        This subroutine performs the iteration described in M. Harrison's
        and Kukushkin's analytic ITER divertor model.
        Report ITER-IL-PH-13-9-e12
        AEA FUS 251: A User's Guide to the PROCESS Systems Code

        :param adas: divertor flux area / main plasma area (long separatrix)
        :type adas: float

        :param aion: ion mass (assumes fuel only) (AMU)
        :type aion: float


        :param anginc: pol. angle of incidence of field line on plate (rad)
        :type anginc: float

        :param c1div: fitting coefficient for plate temperature
        :type c1div: float

        :param c2div: fitting coefficient for plate temperature
        :type c2div: float

        :param c3div: fitting coefficient for heat load
        :type c3div: float

        :param c4div: fitting coefficient for heat load
        :type c4div: float

        :param c5div: fitting coefficient for 'omlarg'
        :type c5div: float

        :param delld: coeff. for power distribution flow into scrapeoff
        :type delld: float

        :param delne: scrapeoff density by main plasma (10**20 m-3)
        :type delne: float

        :param fdfs: gradient ratio (private flux side/other side) in 'omlarg'
        :type fdfs: float

        :param fififi: coeff. used in sheath energy transfer factor calc.
        :type fififi: float

        :param frgd: separatrix area to divertor / total separatrix area
        :type frgd: float

        :param frrp: fraction of radiated power to plate
        :type frrp: float

        :param minstang: minimum strike angle (total) for heat flux calc.
        :type minstang: float

        :param omegan: pressure ratio of (plate / main plasma)
        :type omegan: float

        :param qdiv: heat flux across separatrix to divertor (MW/m2)
        :type qdiv: float

        :param pdiv: power flow to plate (MW)
        :type pdiv: float

        :param rbpbt: ratio of toroidal / poloidal field at strike point
        :type rbpbt: float

        :param rconl: connection length ratio (divertor region/main plasma region)
        :type rconl: float

        :param rmaj: major radius (m)
        :type rmaj: float

        :param rsrd: ratio of separatrix radius / divertor radius
        :type rsrd: float

        :param tconl: connection length along field line by main plasma (m)
        :type tconl: float

        :param xpara: parallel diffusivity in the plasma scrapeoff (m2/s)
        :type xpara: float

        :param xperp: perpend. diffusivity in the plasma scrapeoff (m2/s)
        :type xperp: float

        :returns:
            - delta (`float`) iteration relative error
            - delw (`float`) energy flow thickness in scrape-off (m)
            - dendiv (`float`) plasma density at divertor (10**20 m-3)
            - densin (`float`) peak plasma density at divertor (on separatrix) (10**20 m-3)
            - gamdt (`float`) plasma flow to plate (10**20/s)
            - lamp (`float`) power flow width (m)
            - omlarg (`float`) factor accounting for power flow to private flux region
            - ppdiv (`float`) divertor heat load without radiation (MW/m2)
            - ppdivr (`float`) divertor heat load with radiation (MW/m2)
            - ptpdiv (`float`) peak plasma temperature at the divertor plate (eV)
            - tdiv (`float`) temperature at the plate (eV)
            - tsep (`float`) temperature at the separatrix (eV)
        """

        c27 = 0.2857143e0
        ei = 13.6e0
        epsilon = 0.001e0
        relerr = 1.0e-9

        fprime = c5div * fdfs
        facdenom = fprime * rsrd * (adas / frgd) ** 2 / rconl
        facdenom = max(facdenom, 0.04e0)
        omlarg = 1.0e0 / (rsrd * math.exp(-facdenom))
        omlarg = min(omlarg, 2.0e0)
        coefl = 1.0e0 / delld + rconl / omlarg  # little 'l' in Harrison model

        #  Start iteration on 2 simultaneous equations (Newton's method)

        tdivges = 150.0e0
        tptsges = 0.9e0
        tdiv = tdivges
        tpts = tptsges

        for i in range(15):
            #  Find derivatives for Newton's method

            tptsp = tpts * (1.0e0 + epsilon)
            deltx = tpts * epsilon
            tdivp = tdiv * (1.0e0 + epsilon)
            delty = tdiv * epsilon

            f1 = self.ftpts(
                aion,
                coefl,
                delne,
                fififi,
                omegan,
                omlarg,
                qdiv,
                tconl,
                xpara,
                xperp,
                tpts,
                tdiv,
            )
            f2 = self.ftdiv(
                aion,
                coefl,
                delne,
                fififi,
                omegan,
                omlarg,
                qdiv,
                tconl,
                xpara,
                xperp,
                tpts,
                tdiv,
            )

            f1dx = (
                self.ftpts(
                    aion,
                    coefl,
                    delne,
                    fififi,
                    omegan,
                    omlarg,
                    qdiv,
                    tconl,
                    xpara,
                    xperp,
                    tptsp,
                    tdiv,
                )
                - f1
            ) / deltx
            f1dy = (
                self.ftpts(
                    aion,
                    coefl,
                    delne,
                    fififi,
                    omegan,
                    omlarg,
                    qdiv,
                    tconl,
                    xpara,
                    xperp,
                    tpts,
                    tdivp,
                )
                - f1
            ) / delty
            f2dx = (
                self.ftdiv(
                    aion,
                    coefl,
                    delne,
                    fififi,
                    omegan,
                    omlarg,
                    qdiv,
                    tconl,
                    xpara,
                    xperp,
                    tptsp,
                    tdiv,
                )
                - f2
            ) / deltx
            f2dy = (
                self.ftdiv(
                    aion,
                    coefl,
                    delne,
                    fififi,
                    omegan,
                    omlarg,
                    qdiv,
                    tconl,
                    xpara,
                    xperp,
                    tpts,
                    tdivp,
                )
                - f2
            ) / delty

            denom = f1dx * f2dy - f1dy * f2dx
            if denom == 0.0e0:
                denom = 1.0e-10
            deltpts = (-f2dy * f1 + f1dy * f2) / denom
            deltdiv = (f2dx * f1 - f1dx * f2) / denom

            #  New guess

            tdiv = tdiv + deltdiv
            tpts = tpts + deltpts
            delta = abs(deltdiv / tdiv + deltpts / tpts)

            if delta < relerr:
                break

        tdiv = max(tdiv, 0.1000e0)
        tpts = max(tpts, 0.0010e0)
        tpts = min(tpts, 0.9999e0)

        #  Some other quantities

        ct = max(0.1e0, (c1div + c2div / (tdiv)))
        ptpdiv = tdiv * ct
        gamdiv = self.gammash(fififi, tdiv)  # sheath coefficient
        dendiv = delne / (omegan * tpts)
        eier = self.erprcy(
            tdiv, dendiv
        )  # ionization + radiation energy / recycle event

        tsep = (
            251.0e0
            * (
                (qdiv * tconl) ** 2
                / (c27 * xpara * (1.0e0 - tpts**3.5e0))
                * coefl
                / (xperp * delne)
            )
            ** 0.2222222e0
        )

        cp = max(0.1e0, (c3div + c4div / (tdiv)))
        angle = math.sin(anginc) * rbpbt

        if minstang != 0.0e0:
            angle = max(angle, (minstang / 57.3e0))

        ppdiv = (
            2.48e2
            * (qdiv) ** 1.55556e0
            / (xperp * delne) ** 0.777778e0
            * (c27 * xpara) ** 0.2222222e0
            * tconl**0.555556e0
            * ((1.0e0 - tpts**3.5e0) / coefl) ** 0.222222e0
            / omlarg
            * (1.0e0 + ei / (gamdiv * tdiv))
            / (1.0e0 + eier / (gamdiv * tdiv))
            * angle
            * cp
        )
        ppdivr = ppdiv * (1.0e0 + frrp * (eier - ei) / (gamdiv * tdiv))
        gamdt = 6.25e4 * ppdiv / (gamdiv * ptpdiv)
        densin = omegan * tsep * delne / ptpdiv
        delw = (
            4.01e-3
            * (delne * xperp) ** 0.7777778e0
            * tconl**0.4444444e0
            * coefl**0.2222222e0
            / (
                (qdiv) ** 0.55555556e0
                * (c27 * xpara * (1.0e0 - tpts**3.5e0)) ** 0.22222e0
            )
        )
        lamp = pdiv * rsrd / (2.0e0 * constants.pi * rmaj * ppdiv)

        return (
            delta,
            delw,
            dendiv,
            densin,
            gamdt,
            lamp,
            omlarg,
            ppdiv,
            ppdivr,
            ptpdiv,
            tdiv,
            tsep,
        )

    def erprcy(self, tdiv: float, ndiv: float) -> float:
        """Function providing the (energy radiated + ionized) per neutral
        recycle event from the Harrison / Kukushkin ITER model
        author: J Galambos, ORNL
        author: P J Knight, CCFE, Culham Science Centre


        This function calculates the total energy (in eV) radiated and
        ionized, per neutral recycle event, from the Harrison / Kukushkin
        ITER model.
        Report ITER-IL-PH-13-9-e12
        AEA FUS 251: A User's Guide to the PROCESS Systems Code

        :param tdiv: electron temperature at the plate (eV)
        :type tdiv: float

        :param ndiv: electron density at the plate (10**20 m-3)
        :type ndiv: float

        :returns: the total energy (in eV) radiated and ionized, per neutral recycle event, from the Harrison / Kukushkin ITER model.
        :rtype: float
        """

        return max(17.5e0 + (5.0e0 + 37.5e0 / tdiv) * math.log10(10.0e0 / ndiv), 0.001)

    def ftdiv(
        self,
        aion: float,
        coefl: float,
        delne: float,
        fififi: float,
        omegan: float,
        omlarg: float,
        qdiv: float,
        tconl: float,
        xpara: float,
        xperp: float,
        xx: float,
        yy: float,
    ) -> float:
        """Function for divertor temperature solution
        author: J Galambos, ORNL
        author: P J Knight, CCFE, Culham Science Centre

        This function calculates an estimate for the divertor temperature (eV).
        Report ITER-IL-PH-13-9-e12
        AEA FUS 251: A User's Guide to the PROCESS Systems Code

        :param aion: ion mass (assumes fuel only) (AMU)
        :type aion: float

        :param coefl: little 'l' in Harrison model
        :type coefl: float

        :param delne: scrapeoff density by main plasma (10**20 m-3)
        :type delne: float

        :param fififi: coeff. used in sheath energy transfer factor calc.
        :type fififi: float

        :param omegan: pressure ratio of (plate / main plasma)
        :type omegan: float

        :param omlarg: factor accounting for power flow
        :type omlarg: float

        :param qdiv: heat flux across separatrix to divertor (MW/m2)
        :type qdiv: float

        :param tconl: connection length along field line by main plasma (m)
        :type tconl: float

        :param xpara: parallel diffusivity in the plasma scrapeoff (m2/s)
        :type xpara: float

        :param xperp: perpend. diffusivity in the plasma scrapeoff (m2/s)
        :type xperp: float

        :param xx: T_plate / T_separatrix guess
        :type xx: float

        :param yy: T_plate guess (eV)
        :type yy: float

        :returns: an estimate for the divertor temperature (eV)
        :rtype: float
        """

        c27 = 0.28571428e0

        xxs = max(xx, 0.001e0)
        xxs = min(xxs, 0.99999e0)
        yys = max(yy, 0.1e0)

        dendiv = delne * omegan / xxs
        gamdiv = self.gammash(fififi, yys)
        eier = self.erprcy(yys, dendiv)
        ff = (
            20.16e0
            * aion
            * (
                (qdiv) ** 10
                * (c27 * xpara) ** 4
                / (xperp**5 * delne**14)
                * tconl
                * (1.0e0 - xxs**3.5e0) ** 4
                / coefl**4
            )
            ** 0.22222e0
            / (omegan * gamdiv * omlarg * (1.0e0 + eier / (gamdiv * yys))) ** 2
        )

        return yys - ff

    def ftpts(
        self,
        aion: float,
        coefl: float,
        delne: float,
        fififi: float,
        omegan: float,
        omlarg: float,
        qdiv: float,
        tconl: float,
        xpara: float,
        xperp: float,
        xx: float,
        yy: float,
    ) -> float:
        """Function for divertor model temperature ratio solution
        author: J Galambos, ORNL
        author: P J Knight, CCFE, Culham Science Centre

        This function updates the divertor model temperature ratio solution.
        Report ITER-IL-PH-13-9-e12
        AEA FUS 251: A User's Guide to the PROCESS Systems Code

        :param aion: ion mass (assumes fuel only) (AMU)
        :type aion: float

        :param coefl: little 'l' in Harrison model
        :type coefl: float

        :param delne: scrapeoff density by main plasma (10**20 m-3)
        :type delne: float

        :param fififi: coeff. used in sheath energy transfer factor calc.
        :type fififi: float

        :param omegan: pressure ratio of (plate / main plasma)
        :type omegan: float

        :param omlarg: factor accounting for power flow
        :type omlarg: float

        :param qdiv: heat flux across separatrix to divertor (MW/m2)
        :type qdiv: float

        :param tconl: connection length along field line by main plasma (m)
        :type tconl: float

        :param xpara: parallel diffusivity in the plasma scrapeoff (m2/s)
        :type xpara: float

        :param xperp: perpend. diffusivity in the plasma scrapeoff (m2/s)
        :type xperp: float

        :param xx: T_plate / T_separatrix guess
        :type xx: float

        :param yy: T_plate guess (eV)
        :type yy: float

        :returns: an estimate for the divertor temperature (eV)
        :rtype: float
        """

        xxs = max(xx, 0.001e0)
        xxs = min(xxs, 0.99999e0)
        yys = max(yy, 0.1e0)

        dendiv = delne * omegan / xxs
        gamdiv = self.gammash(fififi, yys)
        eier = self.erprcy(yys, dendiv)

        ff = (
            xxs**3.5e0
            + 9.66e0
            * (xxs / aion) ** 0.9e0
            * (xperp / (qdiv) ** 2) ** 0.8e0
            * coefl
            / (2.0e0 * xpara / 7.0e0)
            * tconl**0.2e0
            * (omlarg * omegan * gamdiv * (1.0e0 + eier / (gamdiv * yys))) ** 1.8e0
            * delne**2.6e0
        )

        return 1.0e0 - ff

    def gammash(self, gcoef: float, tdiv: float) -> float:
        """Function to provide the plasma sheath energy transfer coefficient
        author: J Galambos, ORNL
        author: P J Knight, CCFE, Culham Science Centre

        This function provides the plasma sheath energy transfer coefficient
        from the Harrison / Kukushkin ITER model.
        Report ITER-IL-PH-13-9-e12
        AEA FUS 251: A User's Guide to the PROCESS Systems Code

        :param gcoef: coefficient
        :type gcoef: float

        :param tdiv: electron temperature at the plate (eV)
        :type tdiv: float

        :returns: plasma sheath energy transfer coefficient
        :rtype: float
        """

        return 8.3e0 - 6.0e0 * (0.07e0 - 0.18e0 * math.log10(3.0e0 * tdiv * gcoef))

    def divtart(
        self,
        rmajor: float,
        rminor: float,
        triang: float,
        scrapli: float,
        vgap: float,
        pdivt: float,
        output: bool,
    ) -> float:
        """Tight aspect ratio tokamak divertor model
        author: P J Knight, CCFE, Culham Science Centre

        This subroutine calculates the divertor heat load for a tight aspect
        ratio machine, by assuming that the power is evenly spread around the
        divertor chamber by the action of a gaseous target. Each divertor is
        assumed to be approximately triangular in the R,Z plane.
        AEA FUS 64: Figure 2
        AEA FUS 251: A User's Guide to the PROCESS Systems Code

        :param rmajor: plasma major radius (m)
        :type rmajor: float

        :param rminor: plasma minor radius (m)
        :type rminor: float

        :param triang: plasma triangularity
        :type triang: float

        :param scrapli: inboard scrape-off width (m)
        :type scrapli: float

        :param vgap: top scrape-off width (m)
        :type vgap: float

        :param pdivt: power to the divertor (MW)
        :type pdivt: float

        :returns: divertor heat load for a tight aspect ratio machine
        :rtype: float

        :param output: indicate whether output should be written to the output file, or not
        :type output: boolean
        """

        #  Thickness of centrepost + first wall at divertor height

        r1 = rmajor - rminor * triang - 3.0e0 * scrapli + tfv.drtop

        #  Outer radius of divertor region

        r2 = rmajor + rminor

        #  Angle of diagonal divertor plate from horizontal

        if vgap <= 0.0e0:
            eh.fdiags[0] = vgap
            eh.report_error(22)

        theta = math.atan(vgap / (r2 - r1))

        #  Vertical plate area

        a1 = 2.0e0 * constants.pi * r1 * vgap

        #  Horizontal plate area

        a2 = constants.pi * (r2 * r2 - r1 * r1)

        #  Diagonal plate area

        a3 = a2 / (math.cos(theta) * math.cos(theta))

        #  Total divertor area (N.B. there are two of them)

        areadv = 2.0e0 * (a1 + a2 + a3)
        if dv.i_hldiv == 1:
            dv.hldiv = pdivt / areadv

        if output and dv.i_hldiv == 1:
            po.osubhd(self.outfile, "Divertor Heat Load")
            po.ocmmnt(self.outfile, "Assume an expanded divertor with a gaseous target")
            po.oblnkl(self.outfile)
            po.ovarre(self.outfile, "Power to the divertor (MW)", "(pdivt.)", pdivt)
            po.ovarre(self.outfile, "Divertor surface area (m2)", "(areadv)", areadv)
            po.ovarre(self.outfile, "Divertor heat load (MW/m2)", "(hldiv)", dv.hldiv)

        elif output:
            po.osubhd(self.outfile, "Divertor Heat Load")
            po.ovarre(self.outfile, "Power to the divertor (MW)", "(pdivt.)", pdivt)
            po.ovarre(self.outfile, "Divertor heat load (MW/m2)", "(hldiv)", dv.hldiv)
        return dv.hldiv

    def divwade(
        self,
        rmajor: float,
        rminor: float,
        aspect: float,
        bt: float,
        bp: float,
        pdivt: float,
        flux_exp: float,
        nesep: float,
        beta_div: float,
        rad_fraction_sol: float,
        ftar: float,
        output: bool,
    ) -> float:
        """Divertor heat load model (Wade 2020)
        author: J A Foster, CCFE, Culham Science Centre

        This subroutine calculates the divertor heat flux for any machine,
        with either a single null or double null configuration.
        It uses the Eich scaling (Eich et al. 2013) and spreading factor (Scarabosio et al. 2014)
        to calculate the SOL width. This is then used with a flux expansion factor to calculate
        the wetted area and then the heat flux.

        :param rmajor: plasma major radius (m)
        :type rmajor: float

        :param rminor: plasma minor radius (m)
        :type rminor: float

        :param aspect: tokamak aspect ratio
        :type aspect: float

        :param bt: toroidal field (T)
        :type bt: float

        :param bp: poloidal field (T)
        :type bp: float

        :param pdivt: power to divertor (MW)
        :type pdivt: float

        :param flux_exp: plasma flux expansion in divertor
        :type flux_exp: float

        :param nesep: electron density at separatrix (m-3)
        :type nesep: float

        :param beta_div: field line angle wrt divertor target plate (degrees)
        :type beta_div: float

        :param rad_fraction_sol: SOL radiation fraction
        :type rad_fraction_sol: float

        :param ftar: fraction of power to the lower divertor in double null configuration
        :type ftar: float

        :returns: divertor heat load for a tight aspect ratio machine
        :rtype: float

        :param output: indicate whether output should be written to the output file, or not
        :type output: boolean
        """

        # Radius on midplane
        r_omp = rmajor + rminor

        # B fields on midplane
        Bp_omp = -bp * rmajor / r_omp

        Bt_omp = -bt * rmajor / r_omp

        # Eich scaling for lambda_q
        lambda_eich = (
            1.35 * pdivt**-0.02 * rmajor**0.04 * bp**-0.92 * aspect**0.42
        )

        # Spreading factor
        spread_fact = (
            0.12
            * (nesep / 1e19) ** -0.02
            * pdivt**-0.21
            * rmajor**0.71
            * bp**-0.82
        )

        # SOL width
        lambda_int = lambda_eich + 1.64 * spread_fact

        # Flux angle on midplane
        alpha_mid = math.degrees(math.atan(Bp_omp / Bt_omp))

        # Flux angle in the divertor
        alpha_div = flux_exp * alpha_mid

        # Tilt of the separatrix relative to the target in the poloidal plane
        theta_div = math.asin(
            (1 + 1 / alpha_div**2) * math.sin(math.radians(beta_div))
        )

        # Wetted area
        area_wetted = (
            2 * constants.pi * rmajor * lambda_int * flux_exp * math.sin(theta_div)
        )

        # Divertor heat load
        hldiv_base = pdivt * (1 - rad_fraction_sol) / area_wetted

        # For double null, calculate heat loads to upper and lower divertors and use the highest
        if pv.idivrt == 2:
            hldiv_lower = ftar * hldiv_base
            hldiv_upper = (1.0 - ftar) * hldiv_base
            dv.hldiv = max(hldiv_lower, hldiv_upper)
        else:
            dv.hldiv = hldiv_base

        if output:
            po.osubhd(self.outfile, "Divertor Heat Load")
            po.ocmmnt(self.outfile, "Assume an expanded divertor with a gaseous target")
            po.oblnkl(self.outfile)
            po.ovarre(self.outfile, "Flux expansion", "(flux_exp)", flux_exp)
            po.ovarre(
                self.outfile,
                "Field line angle wrt to target divertor plate (degrees)",
                "(beta_div)",
                beta_div,
            )
            po.ovarre(self.outfile, "Divertor heat load (MW/m2)", "(hldiv)", dv.hldiv)
        return dv.hldiv
