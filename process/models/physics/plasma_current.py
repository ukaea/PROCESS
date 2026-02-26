import logging

import numpy as np

from process import constants
from process.data_structure import (
    physics_variables,
)
from process.exceptions import ProcessValueError

logger = logging.getLogger(__name__)


class PlasmaCurrent:
    """Class to hold plasma current calculations for plasma processing."""

    def __init__(self):
        self.outfile = constants.NOUT
        self.mfile = constants.MFILE

    def calculate_plasma_current(
        self,
        alphaj: float,
        alphap: float,
        b_plasma_toroidal_on_axis: float,
        eps: float,
        i_plasma_current: int,
        kappa: float,
        kappa95: float,
        pres_plasma_on_axis: float,
        len_plasma_poloidal: float,
        q95: float,
        rmajor: float,
        rminor: float,
        triang: float,
        triang95: float,
    ) -> tuple[float, float, float, float, float]:
        """Calculate the plasma current.

        Args:
            alphaj (float): Current profile index.
            alphap (float): Pressure profile index.
            b_plasma_toroidal_on_axis (float): Toroidal field on axis (T).
            eps (float): Inverse aspect ratio.
            i_plasma_current (int): Current scaling model to use.
                1 = Peng analytic fit
                2 = Peng divertor scaling (TART,STAR)
                3 = Simple ITER scaling
                4 = IPDG89 scaling
                5 = Todd empirical scaling I
                6 = Todd empirical scaling II
                7 = Connor-Hastie model
                8 = Sauter scaling (allowing negative triangularity)
                9 = FIESTA ST scaling
            kappa (float): Plasma elongation.
            kappa95 (float): Plasma elongation at 95% surface.
            pres_plasma_on_axis (float): Central plasma pressure (Pa).
            len_plasma_poloidal (float): Plasma perimeter length (m).
            q95 (float): Plasma safety factor at 95% flux (= q-bar for i_plasma_current=2).
            ind_plasma_internal_norm (float): Plasma normalised internal inductance.
            rmajor (float): Major radius (m).
            rminor (float): Minor radius (m).
            triang (float): Plasma triangularity.
            triang95 (float): Plasma triangularity at 95% surface.

        Returns:
            Tuple[float, float, float,]: Tuple containing b_plasma_poloidal_average, qstar, plasma_current,

        Raises:
            ValueError: If invalid value for i_plasma_current is provided.

        Notes:
            This routine calculates the plasma current based on the edge safety factor q95.
            It will also make the current profile parameters consistent with the q-profile if required.

        References:
            - J D Galambos, STAR Code : Spherical Tokamak Analysis and Reactor Code, unpublished internal Oak Ridge document
            - Peng, Y. K. M., Galambos, J. D., & Shipe, P. C. (1992).
              'Small Tokamaks for Fusion Technology Testing'. Fusion Technology, 21(3P2A),
              1729-1738. https://doi.org/10.13182/FST92-A29971
            - ITER Physics Design Guidelines: 1989 [IPDG89], N. A. Uckan et al, ITER Documentation Series No.10, IAEA/ITER/DS/10, IAEA, Vienna, 1990
            - M. Kovari et al, 2014, "PROCESS": A systems code for fusion power plants - Part 1: Physics
            - H. Zohm et al, 2013, On the Physics Guidelines for a Tokamak DEMO
            - T. Hartmann, 2013, Development of a modular systems code to analyse the implications of physics assumptions on the design of a demonstration fusion power plant
            - Sauter, Geometric formulas for systems codes..., FED 2016
        """
        # Aspect ratio
        aspect_ratio = 1.0 / eps

        # Only the Sauter scaling (i_plasma_current=8) is suitable for negative triangularity:
        if i_plasma_current != 8 and triang < 0.0:
            raise ProcessValueError(
                f"Triangularity is negative without i_plasma_current = 8 selected: {triang=}, {i_plasma_current=}"
            )

        # Calculate the function Fq that scales the edge q from the
        # circular cross-section cylindrical case

        # Peng analytical fit
        if i_plasma_current == 1:
            fq = self.current.calculate_current_coefficient_peng(
                eps, len_plasma_poloidal, rminor
            )

        # Peng scaling for double null divertor; TARTs [STAR Code]
        elif i_plasma_current == 2:
            plasma_current = 1.0e6 * self.current.calculate_plasma_current_peng(
                q95, aspect_ratio, eps, rminor, b_plasma_toroidal_on_axis, kappa, triang
            )

        # Simple ITER scaling (simply the cylindrical case)
        elif i_plasma_current == 3:
            fq = 1.0

        # ITER formula (IPDG89)
        elif i_plasma_current == 4:
            fq = self.calculate_current_coefficient_ipdg89(eps, kappa95, triang95)

        # Todd empirical scalings
        # D.C.Robinson and T.N.Todd, Plasma and Contr Fusion 28 (1986) 1181
        elif i_plasma_current in [5, 6]:
            fq = self.calculate_current_coefficient_todd(eps, kappa95, triang95, model=1)

            if i_plasma_current == 6:
                fq = self.calculate_current_coefficient_todd(
                    eps, kappa95, triang95, model=2
                )

        # Connor-Hastie asymptotically-correct expression
        elif i_plasma_current == 7:
            fq = self.calculate_current_coefficient_hastie(
                alphaj,
                alphap,
                b_plasma_toroidal_on_axis,
                triang95,
                eps,
                kappa95,
                pres_plasma_on_axis,
                constants.RMU0,
            )

        # Sauter scaling allowing negative triangularity [FED May 2016]
        # https://doi.org/10.1016/j.fusengdes.2016.04.033.
        elif i_plasma_current == 8:
            # Assumes zero squareness, note takes kappa, delta at separatrix not _95
            fq = self.calculate_current_coefficient_sauter(eps, kappa, triang)

        # FIESTA ST scaling
        # https://doi.org/10.1016/j.fusengdes.2020.111530.
        elif i_plasma_current == 9:
            fq = self.calculate_current_coefficient_fiesta(eps, kappa, triang)
        else:
            raise ProcessValueError(f"Invalid value {i_plasma_current=}")

        # Main plasma current calculation using the fq value from the different settings
        if i_plasma_current != 2:
            plasma_current = (
                (2.0 * np.pi / constants.RMU0)
                * rminor**2
                / (rmajor * q95)
                * fq
                * b_plasma_toroidal_on_axis
            )
        # i_plasma_current == 2 case covered above

        # Calculate cyclindrical safety factor from IPDG89
        qstar = (
            5.0e6
            * rminor**2
            / (rmajor * plasma_current / b_plasma_toroidal_on_axis)
            * 0.5
            * (1.0 + kappa95**2 * (1.0 + 2.0 * triang95**2 - 1.2 * triang95**3))
        )

        # Calculate the poloidal field generated by the plasma current
        b_plasma_poloidal_average = self.calculate_poloidal_field(
            i_plasma_current,
            plasma_current,
            q95,
            aspect_ratio,
            eps,
            b_plasma_toroidal_on_axis,
            kappa,
            triang,
            len_plasma_poloidal,
            constants.RMU0,
        )

        return b_plasma_poloidal_average, qstar, plasma_current

    @staticmethod
    def _plascar_bpol(
        aspect: float, eps: float, kappa: float, delta: float
    ) -> tuple[float, float, float, float]:
        """Calculate the poloidal field coefficients for determining the plasma current
        and poloidal field.


        This internal function calculates the poloidal field coefficients,
        which is used to calculate the poloidal field and the plasma current.

        Parameters
        ----------
        aspect :
            plasma aspect ratio
        eps :
            inverse aspect ratio
        kappa :
            plasma elongation
        delta :
            plasma triangularity

        Returns
        -------
        :
            coefficients ff1, ff2, d1, d2

        References
        ----------
            - Peng, Y. K. M., Galambos, J. D., & Shipe, P. C. (1992).
            'Small Tokamaks for Fusion Technology Testing'. Fusion Technology, 21(3P2A),
            1729-1738. https://doi.org/10.13182/FST92-A29971
            - J D Galambos, STAR Code : Spherical Tokamak Analysis and Reactor Code,
            unpublished internal Oak Ridge document

        """
        # Original coding, only suitable for TARTs [STAR Code]

        c1 = (kappa**2 / (1.0 + delta)) + delta
        c2 = (kappa**2 / (1.0 - delta)) - delta

        d1 = (kappa / (1.0 + delta)) ** 2 + 1.0
        d2 = (kappa / (1.0 - delta)) ** 2 + 1.0

        c1_aspect = ((c1 * eps) - 1.0) if aspect < c1 else (1.0 - (c1 * eps))

        y1 = np.sqrt(c1_aspect / (1.0 + eps)) * ((1.0 + delta) / kappa)
        y2 = np.sqrt((c2 * eps + 1.0) / (1.0 - eps)) * ((1.0 - delta) / kappa)

        h2 = (1.0 + (c2 - 1.0) * (eps / 2.0)) / np.sqrt((1.0 - eps) * (c2 * eps + 1.0))
        f2 = (d2 * (1.0 - delta) * eps) / ((1.0 - eps) * ((c2 * eps) + 1.0))
        g = (eps * kappa) / (1.0 - (eps * delta))
        ff2 = f2 * (g + 2.0 * h2 * np.arctan(y2))

        h1 = (1.0 + (1.0 - c1) * (eps / 2.0)) / np.sqrt((1.0 + eps) * c1_aspect)
        f1 = (d1 * (1.0 + delta) * eps) / ((1.0 + eps) * (c1 * eps - 1.0))

        if aspect < c1:
            ff1 = f1 * (g - h1 * np.log((1.0 + y1) / (1.0 - y1)))
        else:
            ff1 = -f1 * (-g + 2.0 * h1 * np.arctan(y1))

        return ff1, ff2, d1, d2

    def calculate_poloidal_field(
        self,
        i_plasma_current: int,
        ip: float,
        q95: float,
        aspect: float,
        eps: float,
        b_plasma_toroidal_on_axis: float,
        kappa: float,
        delta: float,
        perim: float,
        rmu0: float,
    ) -> float:
        """Function to calculate poloidal field from the plasma current

        This function calculates the poloidal field from the plasma current in Tesla,
        using a simple calculation using Ampere's law for conventional
        tokamaks, or for TARTs, a scaling from Peng, Galambos and
        Shipe (1992).

        Parameters
        ----------
        i_plasma_current :
            current scaling model to use
        ip :
            plasma current (A)
        q95 :
            95% flux surface safety factor
        aspect :
            plasma aspect ratio
        eps :
            inverse aspect ratio
        b_plasma_toroidal_on_axis :
            toroidal field on axis (T)
        kappa :
            plasma elongation
        delta :
            plasma triangularity
        perim :
            plasma perimeter (m)
        rmu0 :
            vacuum permeability (H/m)

        Returns
        -------
        :
            poloidal field in Tesla


        References
        ----------
            - J D Galambos, STAR Code : Spherical Tokamak Analysis and Reactor Code,
            unpublished internal Oak Ridge document
            - Peng, Y. K. M., Galambos, J. D., & Shipe, P. C. (1992).
            'Small Tokamaks for Fusion Technology Testing'. Fusion Technology, 21(3P2A),
            1729-1738. https://doi.org/10.13182/FST92-A29971

        """
        # Use Ampere's law using the plasma poloidal cross-section
        if i_plasma_current != 2:
            return rmu0 * ip / perim
        # Use the relation from Peng, Galambos and Shipe (1992) [STAR code] otherwise
        ff1, ff2, _, _ = self._plascar_bpol(aspect, eps, kappa, delta)

        # Transform q95 to qbar
        qbar = q95 * 1.3e0 * (1.0e0 - physics_variables.eps) ** 0.6e0

        return b_plasma_toroidal_on_axis * (ff1 + ff2) / (2.0 * np.pi * qbar)

    @staticmethod
    def calculate_current_coefficient_peng(
        eps: float, len_plasma_poloidal: float, rminor: float
    ) -> float:
        """
        Calculate the plasma current scaling coefficient for the Peng scaling from the STAR code.

        :param eps: Plasma inverse aspect ratio.
        :type eps: float
        :param len_plasma_poloidal: Plasma poloidal perimeter length [m].
        :type len_plasma_poloidal: float
        :param rminor: Plasma minor radius [m].
        :type rminor: float

        :return: The plasma current scaling coefficient.
        :rtype: float

        :references: None
        """

        return (
            (1.22 - 0.68 * eps)
            / ((1.0 - eps * eps) ** 2)
            * (len_plasma_poloidal / (2.0 * np.pi * rminor)) ** 2
        )

    @staticmethod
    def calculate_plasma_current_peng(
        q95: float,
        aspect: float,
        eps: float,
        rminor: float,
        b_plasma_toroidal_on_axis: float,
        kappa: float,
        delta: float,
    ) -> float:
        """
        Function to calculate plasma current (Peng scaling from the STAR code)

        Parameters:
        - q95: float, 95% flux surface safety factor
        - aspect: float, plasma aspect ratio
        - eps: float, inverse aspect ratio
        - rminor: float, plasma minor radius (m)
        - b_plasma_toroidal_on_axis: float, toroidal field on axis (T)
        - kappa: float, plasma elongation
        - delta: float, plasma triangularity

        Returns:
        - float, plasma current in MA

        This function calculates the plasma current in MA,
        using a scaling from Peng, Galambos and Shipe (1992).
        It is primarily used for Tight Aspect Ratio Tokamaks and is
        selected via i_plasma_current=2.

        References:
        - J D Galambos, STAR Code : Spherical Tokamak Analysis and Reactor Code,
        unpublished internal Oak Ridge document
        - Peng, Y. K. M., Galambos, J. D., & Shipe, P. C. (1992).
        'Small Tokamaks for Fusion Technology Testing'. Fusion Technology, 21(3P2A),
        1729-1738. https://doi.org/10.13182/FST92-A29971
        """

        # Transform q95 to qbar
        qbar = q95 * 1.3e0 * (1.0e0 - physics_variables.eps) ** 0.6e0

        ff1, ff2, d1, d2 = PlasmaCurrent()._plascar_bpol(aspect, eps, kappa, delta)

        e1 = (2.0 * kappa) / (d1 * (1.0 + delta))
        e2 = (2.0 * kappa) / (d2 * (1.0 - delta))

        return (
            rminor
            * b_plasma_toroidal_on_axis
            / qbar
            * 5.0
            * kappa
            / (2.0 * np.pi**2)
            * (np.arcsin(e1) / e1 + np.arcsin(e2) / e2)
            * (ff1 + ff2)
        )

    @staticmethod
    def calculate_current_coefficient_ipdg89(
        eps: float, kappa95: float, triang95: float
    ) -> float:
        """
        Calculate the fq coefficient from the IPDG89 guidlines used in the plasma current scaling.

        Parameters:
        - eps: float, plasma inverse aspect ratio
        - kappa95: float, plasma elongation 95%
        - triang95: float, plasma triangularity 95%

        Returns:
        - float, the fq plasma current coefficient

        This function calculates the fq coefficient used in the IPDG89 plasma current scaling,
        based on the given plasma parameters.

        References:
        - N.A. Uckan and ITER Physics Group, 'ITER Physics Design Guidelines: 1989'
        - T.C.Hender et.al., 'Physics Assesment of the European Reactor Study', AEA FUS 172, 1992
        """
        return (
            0.5
            * (1.17 - 0.65 * eps)
            / ((1.0 - eps * eps) ** 2)
            * (1.0 + kappa95**2 * (1.0 + 2.0 * triang95**2 - 1.2 * triang95**3))
        )

    @staticmethod
    def calculate_current_coefficient_todd(
        eps: float, kappa95: float, triang95: float, model: int
    ) -> float:
        """
        Calculate the fq coefficient used in the two Todd plasma current scalings.

        Parameters:
        - eps: float, plasma inverse aspect ratio
        - kappa95: float, plasma elongation 95%
        - triang95: float, plasma triangularity 95%

        Returns:
        - float, the fq plasma current coefficient

        This function calculates the fq coefficient based on the given plasma parameters for the two Todd scalings.

        References:
        - D.C.Robinson and T.N.Todd, Plasma and Contr Fusion 28 (1986) 1181
        - T.C.Hender et.al., 'Physics Assesment of the European Reactor Study', AEA FUS 172, 1992
        """
        # Calculate the Todd scaling based on the model
        base_scaling = (
            (1.0 + 2.0 * eps**2)
            * ((1.0 + kappa95**2) / 0.5)
            * (
                1.24
                - 0.54 * kappa95
                + 0.3 * (kappa95**2 + triang95**2)
                + 0.125 * triang95
            )
        )
        if model == 1:
            return base_scaling
        if model == 2:
            return base_scaling * (1.0 + (abs(kappa95 - 1.2)) ** 3)
        raise ProcessValueError(f"model = {model} is an invalid option")

    @staticmethod
    def calculate_current_coefficient_hastie(
        alphaj: float,
        alphap: float,
        b_plasma_toroidal_on_axis: float,
        delta95: float,
        eps: float,
        kappa95: float,
        pres_plasma_on_axis: float,
        rmu0: float,
    ) -> float:
        """
        Routine to calculate the f_q coefficient for the Connor-Hastie model used for scaling the plasma current.

        Parameters:
        - alphaj: float, the current profile index
        - alphap: float, the pressure profile index
        - b_plasma_toroidal_on_axis: float, the toroidal field on axis (T)
        - delta95: float, the plasma triangularity 95%
        - eps: float, the inverse aspect ratio
        - kappa95: float, the plasma elongation 95%
        - pres_plasma_on_axis: float, the central plasma pressure (Pa)
        - rmu0: float, the vacuum permeability (H/m)

        Returns:
        - float, the F coefficient

        This routine calculates the f_q coefficient used for scaling the plasma current,
        using the Connor-Hastie scaling

        Reference:
        - J.W.Connor and R.J.Hastie, Culham Lab Report CLM-M106 (1985).
        https://scientific-publications.ukaea.uk/wp-content/uploads/CLM-M106-1.pdf
        - T.C.Hender et.al., 'Physics Assesment of the European Reactor Study', AEA FUS 172, 1992
        """
        # Exponent in Connor-Hastie current profile
        lamda = alphaj

        # Exponent in Connor-Hastie pressure profile
        nu = alphap

        # Central plasma beta
        beta0 = 2.0 * rmu0 * pres_plasma_on_axis / (b_plasma_toroidal_on_axis**2)

        # Plasma internal inductance
        lamp1 = 1.0 + lamda
        li = lamp1 / lamda * (lamp1 / lamda * np.log(lamp1) - 1.0)

        # T/r in AEA FUS 172
        kap1 = kappa95 + 1.0
        tr = kappa95 * delta95 / kap1**2

        # E/r in AEA FUS 172
        er = (kappa95 - 1.0) / kap1

        # T primed in AEA FUS 172
        tprime = 2.0 * tr * lamp1 / (1.0 + 0.5 * lamda)

        # E primed in AEA FUS 172
        eprime = er * lamp1 / (1.0 + lamda / 3.0)

        # Delta primed in AEA FUS 172
        deltap = (0.5 * kap1 * eps * 0.5 * li) + (
            beta0 / (0.5 * kap1 * eps)
        ) * lamp1**2 / (1.0 + nu)

        # Delta/R0 in AEA FUS 172
        deltar = beta0 / 6.0 * (1.0 + 5.0 * lamda / 6.0 + 0.25 * lamda**2) + (
            0.5 * kap1 * eps
        ) ** 2 * 0.125 * (1.0 - (lamda**2) / 3.0)

        # F coefficient
        return (0.5 * kap1) ** 2 * (
            1.0
            + eps**2 * (0.5 * kap1) ** 2
            + 0.5 * deltap**2
            + 2.0 * deltar
            + 0.5 * (eprime**2 + er**2)
            + 0.5 * (tprime**2 + 4.0 * tr**2)
        )

    @staticmethod
    def calculate_current_coefficient_sauter(
        eps: float,
        kappa: float,
        triang: float,
    ) -> float:
        """
        Routine to calculate the f_q coefficient for the Sauter model used for scaling the plasma current.

        Parameters:
        - eps: float, inverse aspect ratio
        - kappa: float, plasma elongation at the separatrix
        - triang: float, plasma triangularity at the separatrix

        Returns:
        - float, the fq coefficient

        Reference:
        - O. Sauter, Geometric formulas for system codes including the effect of negative triangularity,
        Fusion Engineering and Design, Volume 112, 2016, Pages 633-645,
        ISSN 0920-3796, https://doi.org/10.1016/j.fusengdes.2016.04.033.
        """
        w07 = 1.0  # zero squareness - can be modified later if required

        return (
            (4.1e6 / 5.0e6)
            * (1.0 + 1.2 * (kappa - 1.0) + 0.56 * (kappa - 1.0) ** 2)
            * (1.0 + 0.09 * triang + 0.16 * triang**2)
            * (1.0 + 0.45 * triang * eps)
            / (1.0 - 0.74 * eps)
            * (1.0 + 0.55 * (w07 - 1.0))
        )

    @staticmethod
    def calculate_current_coefficient_fiesta(
        eps: float, kappa: float, triang: float
    ) -> float:
        """
        Calculate the fq coefficient used in the FIESTA plasma current scaling.

        Parameters:
        - eps: float, plasma inverse aspect ratio
        - kappa: float, plasma elongation at the separatrix
        - triang: float, plasma triangularity at the separatrix

        Returns:
        - float, the fq plasma current coefficient

        This function calculates the fq coefficient based on the given plasma parameters for the FIESTA scaling.

        References:
        - S.Muldrew et.al,“PROCESS”: Systems studies of spherical tokamaks, Fusion Engineering and Design,
        Volume 154, 2020, 111530, ISSN 0920-3796, https://doi.org/10.1016/j.fusengdes.2020.111530.
        """
        return 0.538 * (1.0 + 2.440 * eps**2.736) * kappa**2.154 * triang**0.060
