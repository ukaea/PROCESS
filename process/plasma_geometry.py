import logging
import numpy
from process.fortran import constants
from process.fortran import build_variables
from process.fortran import physics_variables

logger = logging.getLogger(__name__)


class PlasmaGeom:
    def __init__(self):
        self.outfile = constants.nout

    def geomty(self):
        """
        Plasma geometry parameters
        author: P J Knight, CCFE, Culham Science Centre
        None
        This subroutine calculates the plasma geometry parameters.
        J D Galambos, STAR Code : Spherical Tokamak Analysis and Reactor Code,
        unpublished internal Oak Ridge document
        F/MI/PJK/LOGBOOK14, pp.41-43
        AEA FUS 251: A User's Guide to the PROCESS Systems Code
        H. Zohm et al, On the Physics Guidelines for a Tokamak DEMO,
        FTP/3-3, Proc. IAEA Fusion Energy Conference, October 2012, San Diego
        """

        xsi = 0.0e0
        xso = 0.0e0
        thetai = 0.0e0
        thetao = 0.0e0
        xi = 0.0e0
        xo = 0.0e0

        physics_variables.rminor = physics_variables.rmajor / physics_variables.aspect
        physics_variables.eps = 1.0e0 / physics_variables.aspect

        # Calculate shaping terms, rather than use input values

        if (
            physics_variables.ishape == 0
        ):  # Use input kappa, physics_variables.triang values

            #  Rough estimate of 95% values
            #  ITER Physics Design Guidlines: 1989 (Uckan et al. 1990)
            #  (close to previous estimate of (physics_variables.kappa - 0.04) / 1.1
            #  over a large physics_variables.kappa range)

            physics_variables.kappa95 = physics_variables.kappa / 1.12e0
            physics_variables.triang95 = physics_variables.triang / 1.50e0

        if (
            physics_variables.ishape == 1
        ):  # ST scaling with physics_variables.aspect ratio [STAR Code]

            physics_variables.qlim = 3.0e0 * (
                1.0e0 + 2.6e0 * physics_variables.eps**2.8e0
            )

            physics_variables.kappa = 2.05e0 * (
                1.0e0 + 0.44e0 * physics_variables.eps**2.1e0
            )
            physics_variables.triang = 0.53e0 * (
                1.0e0 + 0.77e0 * physics_variables.eps**3
            )

            # SIM 10/09/2020: Switched to FIESTA ST scaling  from IPDG89
            physics_variables.kappa95 = (
                physics_variables.kappa - 0.39467e0
            ) / 0.90698e0  # Fit to FIESTA (Issue #1086)
            physics_variables.triang95 = (
                physics_variables.triang - 0.048306e0
            ) / 1.3799e0

        if (
            physics_variables.ishape == 2
        ):  # Zohm et al. ITER scaling for elongation, input physics_variables.triang

            physics_variables.kappa = physics_variables.fkzohm * min(
                2.0e0, 1.5e0 + 0.5e0 / (physics_variables.aspect - 1.0e0)
            )

            # ITER Physics Design Guidlines: 1989 (Uckan et al. 1990)
            physics_variables.kappa95 = physics_variables.kappa / 1.12e0
            physics_variables.triang95 = physics_variables.triang / 1.50e0

        if (
            physics_variables.ishape == 3
        ):  # Zohm et al. ITER scaling for elongation, input physics_variables.triang95

            physics_variables.kappa = physics_variables.fkzohm * min(
                2.0e0, 1.5e0 + 0.5e0 / (physics_variables.aspect - 1.0e0)
            )

            # ITER Physics Design Guidlines: 1989 (Uckan et al. 1990)
            physics_variables.triang = 1.5e0 * physics_variables.triang95

            physics_variables.kappa95 = physics_variables.kappa / 1.12e0

        if (
            physics_variables.ishape == 4
        ):  # Use input kappa95, physics_variables.triang95 values

            # ITER Physics Design Guidlines: 1989 (Uckan et al. 1990)
            physics_variables.kappa = 1.12e0 * physics_variables.kappa95
            physics_variables.triang = 1.5e0 * physics_variables.triang95

        if (
            physics_variables.ishape == 5
        ):  # Use input kappa95, physics_variables.triang95 values

            # Fit to MAST data (Issue #1086)
            physics_variables.kappa = 0.91300e0 * physics_variables.kappa95 + 0.38654e0
            physics_variables.triang = (
                0.77394e0 * physics_variables.triang95 + 0.18515e0
            )

        if (
            physics_variables.ishape == 6
        ):  # Use input kappa, physics_variables.triang values

            # Fit to MAST data (Issue #1086)
            physics_variables.kappa95 = (
                physics_variables.kappa - 0.38654e0
            ) / 0.91300e0
            physics_variables.triang95 = (
                physics_variables.triang - 0.18515e0
            ) / 0.77394e0

        if (
            physics_variables.ishape == 7
        ):  # Use input kappa95, physics_variables.triang95 values

            # Fit to FIESTA (Issue #1086)
            physics_variables.kappa = 0.90698e0 * physics_variables.kappa95 + 0.39467e0
            physics_variables.triang = (
                1.3799e0 * physics_variables.triang95 + 0.048306e0
            )

        if (
            physics_variables.ishape == 8
        ):  # Use input kappa, physics_variables.triang values

            # Fit to FIESTA (Issue #1086)
            physics_variables.kappa95 = (
                physics_variables.kappa - 0.39467e0
            ) / 0.90698e0
            physics_variables.triang95 = (
                physics_variables.triang - 0.048306e0
            ) / 1.3799e0

        if (
            physics_variables.ishape == 9
        ):  # Use input triang, physics_variables.rli values

            # physics_variables.kappa found from physics_variables.aspect ratio and plasma internal inductance li(3)
            physics_variables.kappa = (1.09e0 + 0.26e0 / physics_variables.rli) * (
                1.5e0 / physics_variables.aspect
            ) ** 0.4e0

            physics_variables.kappa95 = physics_variables.kappa / 1.12e0
            physics_variables.triang95 = physics_variables.triang / 1.50e0

        if physics_variables.ishape == 10:

            # physics_variables.kappa95 found from physics_variables.aspect ratio and stabilty margin
            # Based on fit to CREATE data. ref Issue #1399
            # valid for EU-DEMO like machine - physics_variables.aspect ratio 2.6 - 3.6
            # Model updated see Issue #1648
            a = 3.68436807e0
            b = -0.27706527e0
            c = 0.87040251e0
            d = -18.83740952e0
            e = -0.27267618e0
            f = 20.5141261e0

            physics_variables.kappa95 = (
                -d
                - c * physics_variables.aspect
                - numpy.sqrt(
                    (c**2.0e0 - 4.0e0 * a * b) * physics_variables.aspect**2.0e0
                    + (2.0e0 * d * c - 4.0e0 * a * e) * physics_variables.aspect
                    + d**2.0e0
                    - 4.0e0 * a * f
                    + 4.0e0 * a * physics_variables.m_s_limit
                )
            ) / (2.0e0 * a)

            if physics_variables.kappa95 > 1.77:
                ratio = 1.77 / physics_variables.kappa95
                corner_fudge = 0.3 * (physics_variables.kappa95 - 1.77) / ratio
                physics_variables.kappa95 = (
                    physics_variables.kappa95 ** (ratio) + corner_fudge
                )

            physics_variables.kappa = 1.12e0 * physics_variables.kappa95
            physics_variables.triang95 = physics_variables.triang / 1.50e0

        if physics_variables.ishape == 11:

            # See Issue #1439
            # physics_variables.triang is an input
            # physics_variables.kappa found from physics_variables.aspect ratio scaling on p32 of Menard:
            # Menard, et al. "Fusion Nuclear Science Facilities
            # and Pilot Plants Based on the Spherical Tokamak."
            # Nucl. Fusion, 2016, 44.

            physics_variables.kappa = 0.95e0 * (
                1.9e0 + 1.9e0 / physics_variables.aspect**1.4e0
            )
            physics_variables.kappa95 = physics_variables.kappa / 1.12e0
            physics_variables.triang95 = physics_variables.triang / 1.50e0

        #  Scrape-off layer thicknesses
        if physics_variables.iscrp == 0:
            build_variables.scraplo = 0.1e0 * physics_variables.rminor
            build_variables.scrapli = 0.1e0 * physics_variables.rminor

        # Find parameters of arcs describing plasma surfaces
        xi, thetai, xo, thetao = self.xparam(
            physics_variables.rminor,
            physics_variables.kappa,
            physics_variables.triang,
        )
        #  Surface area - inboard and outboard.  These are not given by Sauter but
        #  the outboard area is required by DCLL and divertor
        xsi, xso = self.xsurf(
            physics_variables.rmajor,
            physics_variables.rminor,
            xi,
            thetai,
            xo,
            thetao,
        )
        physics_variables.sareao = xso

        # icurr = 8 specifies use of the Sauter geometry as well as plasma current.
        if physics_variables.icurr == 8:
            (
                physics_variables.pperim,
                physics_variables.sf,
                physics_variables.sarea,
                physics_variables.xarea,
                physics_variables.vol,
            ) = self.Sauter_geometry(
                physics_variables.rminor,
                physics_variables.rmajor,
                physics_variables.kappa,
                physics_variables.triang,
            )

        else:

            #  Poloidal perimeter
            physics_variables.pperim = 2.0e0 * (xo * thetao + xi * thetai)
            physics_variables.sf = physics_variables.pperim / (
                2.0e0 * numpy.pi * physics_variables.rminor
            )

            #  Volume
            physics_variables.vol = physics_variables.cvol * self.xvol(
                physics_variables.rmajor,
                physics_variables.rminor,
                xi,
                thetai,
                xo,
                thetao,
            )

            #  Cross-sectional area
            physics_variables.xarea = self.xsecta(xi, thetai, xo, thetao)

            #  Surface area - sum of inboard and outboard.
            physics_variables.sarea = xsi + xso

    def xparam(self, a, kap, tri):
        """
        Routine to find parameters used for calculating geometrical
        properties for double-null plasmas
        author: P J Knight, CCFE, Culham Science Centre
        a      : input real :  plasma minor radius (m)
        kap    : input real :  plasma separatrix elongation
        tri    : input real :  plasma separatrix triangularity
        xi     : output real : radius of arc describing inboard surface (m)
        thetai : output real : half-angle of arc describing inboard surface
        xo     : output real : radius of arc describing outboard surface (m)
        thetao : output real : half-angle of arc describing outboard surface
        This function finds plasma geometrical parameters, using the
        revolution of two intersecting arcs around the device centreline.
        This calculation is appropriate for plasmas with a separatrix.
        F/MI/PJK/LOGBOOK14, p.42
        F/PL/PJK/PROCESS/CODE/047
        AEA FUS 251: A User's Guide to the PROCESS Systems Code
        """
        t = 1.0e0 - tri
        denomi = (kap**2 - t**2) / (2.0e0 * t)
        thetai = numpy.arctan(kap / denomi)
        xi = a * (denomi + 1.0e0 - tri)

        #  Find radius and half-angle of outboard arc

        n = 1.0e0 + tri
        denomo = (kap**2 - n**2) / (2.0e0 * n)
        thetao = numpy.arctan(kap / denomo)
        xo = a * (denomo + 1.0e0 + tri)

        return xi, thetai, xo, thetao

    def surfa(self, a, r, k, d):
        """
        Plasma surface area calculation
        author: P J Knight, CCFE, Culham Science Centre
        a      : input real :  plasma minor radius (m)
        r      : input real :  plasma major radius (m)
        k      : input real :  plasma separatrix elongation
        d      : input real :  plasma separatrix triangularity
        sa     : output real : plasma total surface area (m2)
        so     : output real : plasma outboard surface area (m2)
        This function finds the plasma surface area, using the
        revolution of two intersecting arcs around the device centreline.
        This calculation is appropriate for plasmas with a separatrix.
        It was the original method in PROCESS.
        J D Galambos, STAR Code : Spherical Tokamak Analysis and Reactor Code,
        unpublished internal Oak Ridge document
        """
        radco = a * (1.0e0 + (k**2 + d**2 - 1.0e0) / (2.0e0 * (1.0e0 + d)))
        b = k * a
        thto = numpy.arcsin(b / radco)
        so = 4.0e0 * numpy.pi * radco * ((r + a - radco) * thto + b)

        #  Inboard side

        radci = a * (1.0e0 + (k**2 + d**2 - 1.0e0) / (2.0e0 * (1.0e0 - d)))
        thti = numpy.arcsin(b / radci)
        si = 4.0e0 * numpy.pi * radci * ((r - a + radci) * thti - b)

        sa = so + si

        return sa, so

    def xsurf(self, rmajor, rminor, xi, thetai, xo, thetao):
        """
        Plasma surface area calculation
        author: P J Knight, CCFE, Culham Science Centre
        rmajor : input real :  plasma major radius (m)
        rminor : input real :  plasma minor radius (m)
        xi     : input real :  radius of arc describing inboard surface (m)
        thetai : input real :  half-angle of arc describing inboard surface
        xo     : input real :  radius of arc describing outboard surface (m)
        thetao : input real :  half-angle of arc describing outboard surface
        xsi    : output real : inboard surface area (m2)
        xso    : output real : outboard surface area (m2)
        This function finds the plasma surface area, using the
        revolution of two intersecting arcs around the device centreline.
        This calculation is appropriate for plasmas with a separatrix.
        F/MI/PJK/LOGBOOK14, p.43
        AEA FUS 251: A User's Guide to the PROCESS Systems Code
        """
        fourpi = 4.0e0 * numpy.pi

        rc = rmajor - rminor + xi
        xsi = fourpi * xi * (rc * thetai - xi * numpy.sin(thetai))

        rc = rmajor + rminor - xo
        xso = fourpi * xo * (rc * thetao + xo * numpy.sin(thetao))

        return xsi, xso

    def perim(self, a, kap, tri):
        """
        Plasma poloidal perimeter calculation
        author: P J Knight, CCFE, Culham Science Centre
        a      : input real :  plasma minor radius (m)
        kap    : input real :  plasma separatrix elongation
        tri    : input real :  plasma separatrix triangularity
        This function finds the plasma poloidal perimeter, using the
        revolution of two intersecting arcs around the device centreline.
        This calculation is appropriate for plasmas with a separatrix.
        F/PL/PJK/PROCESS/CODE/047
        AEA FUS 251: A User's Guide to the PROCESS Systems Code
        """

        #  Inboard arc

        denomi = (tri**2 + kap**2 - 1.0e0) / (2.0e0 * (1.0e0 - tri)) + tri
        thetai = numpy.arctan(kap / denomi)
        xli = a * (denomi + 1.0e0 - tri)

        #  Outboard arc

        denomo = (tri**2 + kap**2 - 1.0e0) / (2.0e0 * (1.0e0 + tri)) - tri
        thetao = numpy.arctan(kap / denomo)
        xlo = a * (denomo + 1.0e0 + tri)

        perim = 2.0e0 * (xlo * thetao + xli * thetai)

        return perim

    def xvol(self, rmajor, rminor, xi, thetai, xo, thetao):
        """
        Plasma volume calculation
        author: P J Knight, CCFE, Culham Science Centre
        rmajor : input real :  plasma major radius (m)
        rminor : input real :  plasma minor radius (m)
        xi     : input real :  radius of arc describing inboard surface (m)
        thetai : input real :  half-angle of arc describing inboard surface
        xo     : input real :  radius of arc describing outboard surface (m)
        thetao : input real :  half-angle of arc describing outboard surface
        This function finds the plasma volume, using the
        revolution of two intersecting arcs around the device centreline.
        This calculation is appropriate for plasmas with a separatrix.
        F/MI/PJK/LOGBOOK14, p.43
        AEA FUS 251: A User's Guide to the PROCESS Systems Code
        """
        # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        third = 1.0e0 / 3.0e0

        rc = rmajor - rminor + xi
        vin = (
            constants.twopi
            * xi
            * (
                rc**2 * numpy.sin(thetai)
                - rc * xi * thetai
                - 0.5e0 * rc * xi * numpy.sin(2.0e0 * thetai)
                + xi * xi * numpy.sin(thetai)
                - third * xi * xi * (numpy.sin(thetai)) ** 3
            )
        )

        rc = rmajor + rminor - xo
        vout = (
            constants.twopi
            * xo
            * (
                rc**2 * numpy.sin(thetao)
                + rc * xo * thetao
                + 0.5e0 * rc * xo * numpy.sin(2.0e0 * thetao)
                + xo * xo * numpy.sin(thetao)
                - third * xo * xo * (numpy.sin(thetao)) ** 3
            )
        )

        xvol = vout - vin
        xvol = xvol

        return xvol

    def xsecta(self, xi, thetai, xo, thetao):
        """
        Plasma cross-sectional area calculation
        author: P J Knight, CCFE, Culham Science Centre
        xi     : input real :  radius of arc describing inboard surface (m)
        thetai : input real :  half-angle of arc describing inboard surface
        xo     : input real :  radius of arc describing outboard surface (m)
        thetao : input real :  half-angle of arc describing outboard surface
        This function finds the plasma cross-sectional area, using the
        revolution of two intersecting arcs around the device centreline.
        This calculation is appropriate for plasmas with a separatrix.
        F/MI/PJK/LOGBOOK14, p.41
        AEA FUS 251: A User's Guide to the PROCESS Systems Code
        """

        xsecta = xo**2 * (
            thetao - numpy.cos(thetao) * numpy.sin(thetao)
        ) + xi**2 * (thetai - numpy.cos(thetai) * numpy.sin(thetai))

        return xsecta

    def fvol(self, r, a, kap, tri):
        """
        Plasma volume calculation
        author: P J Knight, CCFE, Culham Science Centre
        r      : input real :  plasma major radius (m)
        a      : input real :  plasma minor radius (m)
        kap    : input real :  plasma separatrix elongation
        tri    : input real :  plasma separatrix triangularity
        This function finds the plasma volume, using the
        revolution of two intersecting arcs around the device centreline.
        This calculation is appropriate for plasmas with a separatrix.
        F/MI/PJK/LOGBOOK14, p.41
        F/PL/PJK/PROCESS/CODE/047
        AEA FUS 251: A User's Guide to the PROCESS Systems Code
        """

        zn = kap * a

        c1 = ((r + a) ** 2 - (r - tri * a) ** 2 - zn**2) / (2.0e0 * (1.0e0 + tri) * a)
        rc1 = r + a - c1
        vout = (
            -0.66666666e0 * numpy.pi * zn**3
            + 2.0e0 * numpy.pi * zn * (c1**2 + rc1**2)
            + 2.0e0
            * numpy.pi
            * c1
            * (zn * numpy.sqrt(rc1**2 - zn**2) + rc1**2 * numpy.arcsin(zn / rc1))
        )

        c2 = (-((r - a) ** 2) + (r - tri * a) ** 2 + zn**2) / (
            2.0e0 * (1.0e0 - tri) * a
        )
        rc2 = c2 - r + a
        vin = (
            -0.66666e0 * numpy.pi * zn**3
            + 2.0e0 * numpy.pi * zn * (rc2**2 + c2**2)
            - 2.0e0
            * numpy.pi
            * c2
            * (zn * numpy.sqrt(rc2**2 - zn**2) + rc2**2 * numpy.arcsin(zn / rc2))
        )

        fvol = vout - vin

        return fvol

    def xsect0(self, a, kap, tri):
        """
        Plasma cross-sectional area calculation
        author: P J Knight, CCFE, Culham Science Centre
        a      : input real :  plasma minor radius (m)
        kap    : input real :  plasma separatrix elongation
        tri    : input real :  plasma separatrix triangularity
        This function finds the plasma cross-sectional area, using the
        revolution of two intersecting arcs around the device centreline.
        This calculation is appropriate for plasmas with a separatrix.
        The method for finding the arc radii and angles are copied from
        routine <A HREF="perim.html">PERIM</A>, and are thought to be
        by Peng.
        F/MI/PJK/LOGBOOK14, p.41
        F/PL/PJK/PROCESS/CODE/047
        AEA FUS 251: A User's Guide to the PROCESS Systems Code
        """

        denomi = (tri**2 + kap**2 - 1.0e0) / (2.0e0 * (1.0e0 - tri)) + tri
        thetai = numpy.arctan(kap / denomi)
        xli = a * (denomi + 1.0e0 - tri)

        cti = numpy.cos(thetai)
        sti = numpy.sin(thetai)

        #  Find radius and half-angle of outboard arc

        denomo = (tri**2 + kap**2 - 1.0e0) / (2.0e0 * (1.0e0 + tri)) - tri
        thetao = numpy.arctan(kap / denomo)
        xlo = a * (denomo + 1.0e0 + tri)

        cto = numpy.cos(thetao)
        sto = numpy.sin(thetao)

        #  Find cross-sectional area

        xsect0 = xlo**2 * (thetao - cto * sto) + xli**2 * (thetai - cti * sti)

        return xsect0

    def sauter_geometry(self, a, r0, kap, tri):
        """
        Plasma geometry based on equations (36) in O. Sauter, Fusion Engineering and Design 112 (2016) 633–645
        'Geometric formulas for system codes including the effect of negative triangularity'
        Author: Michael Kovari, issue #392
        a      : input real :  plasma minor radius (m)
        r0     : input real :  plasma major radius (m)
        kap    : input real :  plasma separatrix elongation
        tri    : input real :  plasma separatrix triangularity
        """
        w07 = 1
        eps = a / r0

        # Poloidal perimeter (named Lp in Sauter)
        pperim = (
            2.0e0
            * numpy.pi
            * a
            * (1 + 0.55 * (kap - 1))
            * (1 + 0.08 * tri**2)
            * (1 + 0.2 * (w07 - 1))
        )

        # A geometric factor
        sf = pperim / (2.0e0 * numpy.pi * a)

        # Surface area (named Ap in Sauter)
        sarea = 2.0e0 * numpy.pi * r0 * (1 - 0.32 * tri * eps) * pperim

        # Cross-section area (named S_phi in Sauter)
        xarea = numpy.pi * a**2 * kap * (1 + 0.52 * (w07 - 1))

        # Volume
        vol = 2.0e0 * numpy.pi * r0 * (1 - 0.25 * tri * eps) * xarea

        return pperim, sf, sarea, xarea, vol
