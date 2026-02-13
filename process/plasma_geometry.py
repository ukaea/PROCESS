import logging

import numpy as np

from process import constants
from process.data_structure import build_variables, physics_variables

logger = logging.getLogger(__name__)


class PlasmaGeom:
    def __init__(self):
        self.outfile = constants.NOUT

    def plasma_geometry(self):
        """Plasma geometry parameters

        This method calculates the plasma geometry parameters based on various shaping terms and input values.
        It updates the `physics_variables` with calculated values for kappa, triangularity, surface area, volume, etc.

        :references:
            - J D Galambos, STAR Code : Spherical Tokamak Analysis and Reactor Code,
              unpublished internal Oak Ridge document
            - H. Zohm et al, On the Physics Guidelines for a Tokamak DEMO,
              FTP/3-3, Proc. IAEA Fusion Energy Conference, October 2012, San Diego
        """

        xsi = 0.0e0
        xso = 0.0e0
        thetai = 0.0e0
        thetao = 0.0e0
        xi = 0.0e0
        xo = 0.0e0

        # Define plasma minor radius from major radius and aspect ratio
        physics_variables.rminor = physics_variables.rmajor / physics_variables.aspect

        # Define the inverse aspect ratio
        physics_variables.eps = 1.0e0 / physics_variables.aspect

        # ======================================================================

        if (
            physics_variables.i_plasma_geometry == 0
        ):  # Use input kappa, physics_variables.triang values
            #  Rough estimate of 95% values
            #  ITER Physics Design Guidlines: 1989 (Uckan et al. 1990)
            #  (close to previous estimate of (physics_variables.kappa - 0.04) / 1.1
            #  over a large physics_variables.kappa range)

            physics_variables.kappa95 = physics_variables.kappa / 1.12e0
            physics_variables.triang95 = physics_variables.triang / 1.50e0

        # ======================================================================

        if (
            physics_variables.i_plasma_geometry == 1
        ):  # ST scaling with physics_variables.aspect ratio [STAR Code]
            physics_variables.q95_min = 3.0e0 * (
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

        # ======================================================================

        if (
            physics_variables.i_plasma_geometry == 2
        ):  # Zohm et al. ITER scaling for elongation, input physics_variables.triang
            physics_variables.kappa = physics_variables.fkzohm * min(
                2.0e0, 1.5e0 + 0.5e0 / (physics_variables.aspect - 1.0e0)
            )

            # ITER Physics Design Guidlines: 1989 (Uckan et al. 1990)
            physics_variables.kappa95 = physics_variables.kappa / 1.12e0
            physics_variables.triang95 = physics_variables.triang / 1.50e0

        # ======================================================================

        if (
            physics_variables.i_plasma_geometry == 3
        ):  # Zohm et al. ITER scaling for elongation, input physics_variables.triang95
            physics_variables.kappa = physics_variables.fkzohm * min(
                2.0e0, 1.5e0 + 0.5e0 / (physics_variables.aspect - 1.0e0)
            )

            # ITER Physics Design Guidlines: 1989 (Uckan et al. 1990)
            physics_variables.triang = 1.5e0 * physics_variables.triang95

            physics_variables.kappa95 = physics_variables.kappa / 1.12e0

        # ======================================================================

        if (
            physics_variables.i_plasma_geometry == 4
        ):  # Use input kappa95, physics_variables.triang95 values
            # ITER Physics Design Guidlines: 1989 (Uckan et al. 1990)
            physics_variables.kappa = 1.12e0 * physics_variables.kappa95
            physics_variables.triang = 1.5e0 * physics_variables.triang95

        # ======================================================================

        if (
            physics_variables.i_plasma_geometry == 5
        ):  # Use input kappa95, physics_variables.triang95 values
            # Fit to MAST data (Issue #1086)
            physics_variables.kappa = 0.91300e0 * physics_variables.kappa95 + 0.38654e0
            physics_variables.triang = 0.77394e0 * physics_variables.triang95 + 0.18515e0

        # ======================================================================

        if (
            physics_variables.i_plasma_geometry == 6
        ):  # Use input kappa, physics_variables.triang values
            # Fit to MAST data (Issue #1086)
            physics_variables.kappa95 = (physics_variables.kappa - 0.38654e0) / 0.91300e0
            physics_variables.triang95 = (
                physics_variables.triang - 0.18515e0
            ) / 0.77394e0

        # ======================================================================

        if (
            physics_variables.i_plasma_geometry == 7
        ):  # Use input kappa95, physics_variables.triang95 values
            # Fit to FIESTA (Issue #1086)
            physics_variables.kappa = 0.90698e0 * physics_variables.kappa95 + 0.39467e0
            physics_variables.triang = 1.3799e0 * physics_variables.triang95 + 0.048306e0

        # ======================================================================

        if (
            physics_variables.i_plasma_geometry == 8
        ):  # Use input kappa, physics_variables.triang values
            # Fit to FIESTA (Issue #1086)
            physics_variables.kappa95 = (physics_variables.kappa - 0.39467e0) / 0.90698e0
            physics_variables.triang95 = (
                physics_variables.triang - 0.048306e0
            ) / 1.3799e0

        # ======================================================================

        if (
            physics_variables.i_plasma_geometry == 9
        ):  # Use input triang, physics_variables.ind_plasma_internal_norm values
            # physics_variables.kappa found from physics_variables.aspect ratio and plasma internal inductance li(3)
            physics_variables.kappa = (
                1.09e0 + 0.26e0 / physics_variables.ind_plasma_internal_norm
            ) * (1.5e0 / physics_variables.aspect) ** 0.4e0

            physics_variables.kappa95 = physics_variables.kappa / 1.12e0
            physics_variables.triang95 = physics_variables.triang / 1.50e0

        # ======================================================================

        if physics_variables.i_plasma_geometry == 10:
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
                - np.sqrt(
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

        # ======================================================================

        if physics_variables.i_plasma_geometry == 11:
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

        # ======================================================================

        #  Scrape-off layer thicknesses
        if physics_variables.i_plasma_wall_gap == 0:
            build_variables.dr_fw_plasma_gap_outboard = 0.1e0 * physics_variables.rminor
            build_variables.dr_fw_plasma_gap_inboard = 0.1e0 * physics_variables.rminor

        # ======================================================================

        # Find parameters of arcs describing plasma surfaces
        xi, thetai, xo, thetao = self.plasma_angles_arcs(
            physics_variables.rminor,
            physics_variables.kappa,
            physics_variables.triang,
        )

        #  Surface area - inboard and outboard.  These are not given by Sauter but
        #  the outboard area is required by DCLL and divertor
        xsi, xso = self.plasma_surface_area(
            physics_variables.rmajor,
            physics_variables.rminor,
            xi,
            thetai,
            xo,
            thetao,
        )
        physics_variables.a_plasma_surface_outboard = xso

        # ======================================================================

        # i_plasma_current = 8 specifies use of the Sauter geometry as well as plasma current.
        if (
            physics_variables.i_plasma_current == 8
            or physics_variables.i_plasma_shape == 1
        ):
            (
                physics_variables.len_plasma_poloidal,
                physics_variables.a_plasma_surface,
                physics_variables.a_plasma_poloidal,
                physics_variables.vol_plasma,
            ) = self.sauter_geometry(
                physics_variables.rminor,
                physics_variables.rmajor,
                physics_variables.kappa,
                physics_variables.triang,
                physics_variables.plasma_square,
            )

        else:
            #  Poloidal perimeter
            physics_variables.len_plasma_poloidal = self.plasma_poloidal_perimeter(
                xi, thetai, xo, thetao
            )

            #  Volume
            physics_variables.vol_plasma = (
                physics_variables.f_vol_plasma
                * self.plasma_volume(
                    physics_variables.rmajor,
                    physics_variables.rminor,
                    xi,
                    thetai,
                    xo,
                    thetao,
                )
            )

            #  Cross-sectional area
            physics_variables.a_plasma_poloidal = self.plasma_cross_section(
                xi, thetai, xo, thetao
            )

            #  Surface area - sum of inboard and outboard.
            physics_variables.a_plasma_surface = xsi + xso

        # ======================================================================

    @staticmethod
    def plasma_angles_arcs(
        a: float, kappa: float, triang: float
    ) -> tuple[float, float, float, float]:
        """Routine to find parameters used for calculating geometrical
        properties for double-null plasmas.

        This function finds plasma geometrical parameters, using the
        revolution of two intersecting arcs around the device centreline.
        This calculation is appropriate for plasmas with a separatrix.


        Parameters
        ----------
        a : float
            Plasma minor radius (m)
        kappa : float
            Plasma separatrix elongation
        tri : float
            Plasma separatrix triangularity

        Returns
        -------
        tuple
            A tuple containing:
            - xi (float): Radius of arc describing inboard surface (m)
            - thetai (float): Half-angle of arc describing inboard surface
            - xo (float): Radius of arc describing outboard surface (m)
            - thetao (float): Half-angle of arc describing outboard surface


        :references:
        - F/MI/PJK/LOGBOOK14, p.42
        - F/PL/PJK/PROCESS/CODE/047
        """
        t = 1.0e0 - triang
        denomi = (kappa**2 - t**2) / (2.0e0 * t)
        thetai = np.arctan(kappa / denomi)
        xi = a * (denomi + 1.0e0 - triang)

        #  Find radius and half-angle of outboard arc

        n = 1.0e0 + triang
        denomo = (kappa**2 - n**2) / (2.0e0 * n)
        thetao = np.arctan(kappa / denomo)
        xo = a * (denomo + 1.0e0 + triang)

        return xi, thetai, xo, thetao

    @staticmethod
    def plasma_poloidal_perimeter(
        xi: float, thetai: float, xo: float, thetao: float
    ) -> float:
        """Calculate the poloidal perimeter of the plasma.

        Parameters
        ----------
        xi : float
            Radius of arc describing inboard surface (m)
        thetai : float
            Half-angle of arc describing inboard surface
        xo : float
            Radius of arc describing outboard surface (m)
        thetao : float
            Half-angle of arc describing outboard surface

        Returns
        -------
        float
            Poloidal perimeter (m)
        """
        return 2.0e0 * (xo * thetao + xi * thetai)

    @staticmethod
    def plasma_surface_area(
        rmajor: float, rminor: float, xi: float, thetai: float, xo: float, thetao: float
    ) -> tuple[float, float]:
        """Plasma surface area calculation

        This function finds the plasma surface area, using the
        revolution of two intersecting arcs around the device centreline.
        This calculation is appropriate for plasmas with a separatrix.

        Parameters
        ----------
        rmajor : float
            Plasma major radius (m)
        rminor : float
            Plasma minor radius (m)
        xi : float
            Radius of arc describing inboard surface (m)
        thetai : float
            Half-angle of arc describing inboard surface
        xo : float
            Radius of arc describing outboard surface (m)
        thetao : float
            Half-angle of arc describing outboard surface

        Returns
        -------
        tuple
            A tuple containing:
            - xsi (float): Inboard surface area (m^2)
            - xso (float): Outboard surface area (m^2)


        :references:
        - F/MI/PJK/LOGBOOK14, p.43
        """
        fourpi = 4.0e0 * np.pi

        rc = rmajor - rminor + xi
        xsi = fourpi * xi * (rc * thetai - xi * np.sin(thetai))

        rc = rmajor + rminor - xo
        xso = fourpi * xo * (rc * thetao + xo * np.sin(thetao))

        return xsi, xso

    @staticmethod
    def plasma_volume(
        rmajor: float, rminor: float, xi: float, thetai: float, xo: float, thetao: float
    ) -> float:
        """Plasma volume calculation
        This function finds the plasma volume, using the
        revolution of two intersecting arcs around the device centreline.
        This calculation is appropriate for plasmas with a separatrix.

        Parameters
        ----------
        rmajor : float
            Plasma major radius (m)
        rminor : float
            Plasma minor radius (m)
        xi : float
            Radius of arc describing inboard surface (m)
        thetai : float
            Half-angle of arc describing inboard surface
        xo : float
            Radius of arc describing outboard surface (m)
        thetao : float
            Half-angle of arc describing outboard surface

        Returns
        -------
        float
            Plasma volume (m^3)


        :references:
        - F/MI/PJK/LOGBOOK14, p.43
        """

        third = 1.0e0 / 3.0e0

        rc = rmajor - rminor + xi
        vin = (
            2.0
            * np.pi
            * xi
            * (
                rc**2 * np.sin(thetai)
                - rc * xi * thetai
                - 0.5e0 * rc * xi * np.sin(2.0e0 * thetai)
                + xi * xi * np.sin(thetai)
                - third * xi * xi * (np.sin(thetai)) ** 3
            )
        )

        rc = rmajor + rminor - xo
        vout = (
            2.0
            * np.pi
            * xo
            * (
                rc**2 * np.sin(thetao)
                + rc * xo * thetao
                + 0.5e0 * rc * xo * np.sin(2.0e0 * thetao)
                + xo * xo * np.sin(thetao)
                - third * xo * xo * (np.sin(thetao)) ** 3
            )
        )

        return vout - vin

    @staticmethod
    def plasma_cross_section(
        xi: float, thetai: float, xo: float, thetao: float
    ) -> float:
        """Plasma cross-sectional area calculation

        This function finds the plasma cross-sectional area, using the
        revolution of two intersecting arcs around the device centreline.
        This calculation is appropriate for plasmas with a separatrix.

        Parameters
        ----------
        xi : float
            Radius of arc describing inboard surface (m)
        thetai : float
            Half-angle of arc describing inboard surface
        xo : float
            Radius of arc describing outboard surface (m)
        thetao : float
            Half-angle of arc describing outboard surface

        Returns
        -------
        float
            Plasma cross-sectional area (m^2)


        :references:
        - F/MI/PJK/LOGBOOK14, p.41
        """

        return xo**2 * (thetao - np.cos(thetao) * np.sin(thetao)) + xi**2 * (
            thetai - np.cos(thetai) * np.sin(thetai)
        )

    @staticmethod
    def sauter_geometry(
        a: float, r0: float, kappa: float, triang: float, square: float
    ) -> tuple[float, float, float, float, float]:
        """Calculate the plasma geometry parameters using the Sauter geometry model.

        Parameters
        ----------
        a : float
            Plasma minor radius (m)
        r0 : float
            Plasma major radius (m)
        kappa : float
            Plasma separatrix elongation
        triang : float
            Plasma separatrix triangularity
        square : float
            Plasma squareness
        a:

        r0:

        kappa:

        triang:

        square:


        Returns
        -------
        tuple
            A tuple containing:
            - len_plasma_poloidal (float): Poloidal perimeter
            - a_plasma_surface (float): Surface area
            - a_plasma_poloidal (float): Cross-section area
            - vol_plasma (float): Plasma volume


        :references:
            - O. Sauter, “Geometric formulas for system codes including the effect of negative triangularity,”
            Fusion Engineering and Design, vol. 112, pp. 633-645, Nov. 2016,
            doi: https://doi.org/10.1016/j.fusengdes.2016.04.033.
        """

        # Calculate w07 parameter from paper from squareness assuming top-down symmetry
        w07 = square + 1

        # Inverse aspect ratio
        eps = a / r0

        # Poloidal perimeter (named Lp in Sauter)
        len_plasma_poloidal = (
            2.0e0
            * np.pi
            * a
            * (1 + 0.55 * (kappa - 1))
            * (1 + 0.08 * triang**2)
            * (1 + 0.2 * (w07 - 1))
        )

        # Surface area (named Ap in Sauter)
        a_plasma_surface = (
            2.0e0 * np.pi * r0 * (1 - 0.32 * triang * eps) * len_plasma_poloidal
        )

        # Cross-section area (named S_phi in Sauter)
        a_plasma_poloidal = np.pi * a**2 * kappa * (1 + 0.52 * (w07 - 1))

        # Volume
        vol_plasma = 2.0e0 * np.pi * r0 * (1 - 0.25 * triang * eps) * a_plasma_poloidal

        return (
            len_plasma_poloidal,
            a_plasma_surface,
            a_plasma_poloidal,
            vol_plasma,
        )


# --------------------------------
# Obsolete legacy calculations
# --------------------------------


def surfa(a: float, r: float, k: float, d: float) -> tuple[float, float]:
    """Plasma surface area calculation

    This function finds the plasma surface area, using the
    revolution of two intersecting arcs around the device centreline.
    This calculation is appropriate for plasmas with a separatrix.
    It was the original method in PROCESS.

    Parameters
    ----------
    a : float
        Plasma minor radius (m)
    r : float
        Plasma major radius (m)
    k : float
        Plasma separatrix elongation
    d : float
        Plasma separatrix triangularity

    Returns
    -------
    tuple
        A tuple containing:
        - sa (float): Plasma total surface area (m^2)
        - so (float): Plasma outboard surface area (m^2)


    :references:
    - J D Galambos, STAR Code : Spherical Tokamak Analysis and Reactor Code,
    unpublished internal Oak Ridge document
    """
    radco = a * (1.0e0 + (k**2 + d**2 - 1.0e0) / (2.0e0 * (1.0e0 + d)))
    b = k * a
    thto = np.arcsin(b / radco)
    so = 4.0e0 * np.pi * radco * ((r + a - radco) * thto + b)

    #  Inboard side

    radci = a * (1.0e0 + (k**2 + d**2 - 1.0e0) / (2.0e0 * (1.0e0 - d)))
    thti = np.arcsin(b / radci)
    si = 4.0e0 * np.pi * radci * ((r - a + radci) * thti - b)

    sa = so + si

    return sa, so


def perim(a: float, kap: float, tri: float) -> float:
    """Plasma poloidal perimeter calculation

    This function finds the plasma poloidal perimeter, using the
    revolution of two intersecting arcs around the device centreline.
    This calculation is appropriate for plasmas with a separatrix.

    Parameters
    ----------
    a : float
        Plasma minor radius (m)
    kap : float
        Plasma separatrix elongation
    tri : float
        Plasma separatrix triangularity
    a:

    kap:

    tri:


    Returns
    -------
    float
        Plasma poloidal perimeter (m)


    :references:
    - F/PL/PJK/PROCESS/CODE/047
    """

    #  Inboard arc

    denomi = (tri**2 + kap**2 - 1.0e0) / (2.0e0 * (1.0e0 - tri)) + tri
    thetai = np.arctan(kap / denomi)
    xli = a * (denomi + 1.0e0 - tri)

    #  Outboard arc

    denomo = (tri**2 + kap**2 - 1.0e0) / (2.0e0 * (1.0e0 + tri)) - tri
    thetao = np.arctan(kap / denomo)
    xlo = a * (denomo + 1.0e0 + tri)

    return 2.0e0 * (xlo * thetao + xli * thetai)


def fvol(r: float, a: float, kap: float, tri: float) -> float:
    """Plasma volume calculation

    This function finds the plasma volume, using the
    revolution of two intersecting arcs around the device centreline.
    This calculation is appropriate for plasmas with a separatrix.

    Parameters
    ----------
    r : float
        Plasma major radius (m)
    a : float
        Plasma minor radius (m)
    kap : float
        Plasma separatrix elongation
    tri : float
        Plasma separatrix triangularity

    Returns
    -------
    float
        Plasma volume (m^3)


    :references:
    - F/MI/PJK/LOGBOOK14, p.41
    - F/PL/PJK/PROCESS/CODE/047
    """

    zn = kap * a

    c1 = ((r + a) ** 2 - (r - tri * a) ** 2 - zn**2) / (2.0e0 * (1.0e0 + tri) * a)
    rc1 = r + a - c1
    vout = (
        -0.66666666e0 * np.pi * zn**3
        + 2.0e0 * np.pi * zn * (c1**2 + rc1**2)
        + 2.0e0
        * np.pi
        * c1
        * (zn * np.sqrt(rc1**2 - zn**2) + rc1**2 * np.arcsin(zn / rc1))
    )

    c2 = (-((r - a) ** 2) + (r - tri * a) ** 2 + zn**2) / (2.0e0 * (1.0e0 - tri) * a)
    rc2 = c2 - r + a
    vin = (
        -0.66666e0 * np.pi * zn**3
        + 2.0e0 * np.pi * zn * (rc2**2 + c2**2)
        - 2.0e0
        * np.pi
        * c2
        * (zn * np.sqrt(rc2**2 - zn**2) + rc2**2 * np.arcsin(zn / rc2))
    )

    return vout - vin


def xsect0(a: float, kap: float, tri: float) -> float:
    """Plasma cross-sectional area calculation

    This function finds the plasma cross-sectional area, using the
    revolution of two intersecting arcs around the device centreline.
    This calculation is appropriate for plasmas with a separatrix.

    Parameters
    ----------
    a : float
        Plasma minor radius (m)
    kap : float
        Plasma separatrix elongation
    tri : float
        Plasma separatrix triangularity

    Returns
    -------
    float
        Plasma cross-sectional area (m^2)

    :references:
    - F/MI/PJK/LOGBOOK14, p.41
    - F/PL/PJK/PROCESS/CODE/047
    """

    denomi = (tri**2 + kap**2 - 1.0e0) / (2.0e0 * (1.0e0 - tri)) + tri
    thetai = np.arctan(kap / denomi)
    xli = a * (denomi + 1.0e0 - tri)

    cti = np.cos(thetai)
    sti = np.sin(thetai)

    #  Find radius and half-angle of outboard arc

    denomo = (tri**2 + kap**2 - 1.0e0) / (2.0e0 * (1.0e0 + tri)) - tri
    thetao = np.arctan(kap / denomo)
    xlo = a * (denomo + 1.0e0 + tri)

    cto = np.cos(thetao)
    sto = np.sin(thetao)

    #  Find cross-sectional area

    return xlo**2 * (thetao - cto * sto) + xli**2 * (thetai - cti * sti)
