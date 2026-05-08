"""Module for plasma L-H and L-I transition power threshold calculations."""

import logging
from enum import IntEnum

from process.core import constants
from process.core import process_output as po
from process.data_structure import (
    numerics,
    physics_variables,
)

logger = logging.getLogger(__name__)


class PlasmaConfinementTransitionModel(IntEnum):
    """Enum for plasma L -> H and L -> I transition power threshold models."""

    ITER1996_NOMINAL = (1, "ITER-1996 Nominal")
    ITER1996_UPPER = (2, "ITER-1996 Upper")
    ITER1996_LOWER = (3, "ITER-1996 Lower")
    SNIPES1997_ITER = (4, "Snipes 1997 ITER Scaling I")
    SNIPES1997_KAPPA = (5, "Snipes 1997 ITER Scaling II")
    MARTIN08_NOMINAL = (6, "Martin 2008 Nominal")
    MARTIN08_UPPER = (7, "Martin 2008 Upper")
    MARTIN08_LOWER = (8, "Martin 2008 Lower")
    SNIPES2000_NOMINAL = (9, "Snipes 2000 Nominal")
    SNIPES2000_UPPER = (10, "Snipes 2000 Upper")
    SNIPES2000_LOWER = (11, "Snipes 2000 Lower")
    SNIPES2000_CLOSED_DIVERTOR_NOMINAL = (12, "Snipes 2000 Closed Divertor Nominal")
    SNIPES2000_CLOSED_DIVERTOR_UPPER = (13, "Snipes 2000 Closed Divertor Upper")
    SNIPES2000_CLOSED_DIVERTOR_LOWER = (14, "Snipes 2000 Closed Divertor Lower")
    HUBBARD2012_NOMINAL = (15, "Hubbard 2012 Nominal")
    HUBBARD2012_LOWER = (16, "Hubbard 2012 Lower")
    HUBBARD2012_UPPER = (17, "Hubbard 2012 Upper")
    HUBBARD2017_I_MODE = (18, "Hubbard 2017 I-Mode")
    MARTIN08_ASPECT_NOMINAL = (19, "Martin 2008 Aspect Corrected Nominal")
    MARTIN08_ASPECT_UPPER = (20, "Martin 2008 Aspect Corrected Upper")
    MARTIN08_ASPECT_LOWER = (21, "Martin 2008 Aspect Corrected Lower")

    def __new__(cls, value: int, full_name: str):
        """Create a new PlasmaConfinementTransitionModel instance.

        Parameters
        ----------
        value : int
            The integer value of the enum member.
        full_name : str
            The full descriptive name of the enum member.

        Returns
        -------
        PlasmaConfinementTransitionModel
            A new instance of PlasmaConfinementTransitionModel.
        """
        obj = int.__new__(cls, value)
        obj._value_ = value
        obj.full_name = full_name
        return obj


class PlasmaConfinementTransition:
    """Class to calculate plasma L -> H and L -> I transition power thresholds."""

    def __init__(self):
        self.outfile = constants.NOUT
        self.mfile = constants.MFILE

    def run(self) -> None:
        """Calculate L- to H-mode and L- to I-mode power thresholds.

        Calculates the plasma L- to H-mode power threshold for different scalings
        and sets the enforced L-H power threshold value based on the selected scaling.
        """
        # Calculate L- to H-mode power threshold for different scalings
        physics_variables.l_h_threshold_powers = self.l_h_threshold_power(
            physics_variables.nd_plasma_electron_line,
            physics_variables.b_plasma_toroidal_on_axis,
            physics_variables.rmajor,
            physics_variables.rminor,
            physics_variables.kappa,
            physics_variables.a_plasma_surface,
            physics_variables.m_ions_total_amu,
            physics_variables.aspect,
            physics_variables.plasma_current,
        )

        # Enforced L-H power threshold value (if constraint 15 is turned on)
        physics_variables.p_l_h_threshold_mw = physics_variables.l_h_threshold_powers[
            physics_variables.i_l_h_threshold - 1
        ]

    def l_h_threshold_power(
        self,
        nd_plasma_electron_line: float,
        b_plasma_toroidal_on_axis: float,
        rmajor: float,
        rminor: float,
        kappa: float,
        a_plasma_surface: float,
        m_ions_total_amu: float,
        aspect: float,
        plasma_current: float,
    ) -> list[float]:
        """L-mode to H-mode power threshold calculation.

        Parameters
        ----------
        nd_plasma_electron_line : float
            Line-averaged electron density (/m3)
        b_plasma_toroidal_on_axis : float
            Toroidal field on axis (T)
        rmajor : float
            Plasma major radius (m)
        rminor : float
            Plasma minor radius (m)
        kappa : float
            Plasma elongation
        a_plasma_surface : float
            Plasma surface area (m2)
        m_ions_total_amu : float
            Average mass of all ions (amu)
        aspect : float
            Aspect ratio
        plasma_current : float
            Plasma current (A)

        Returns
        -------
        list[float]
            Array of power thresholds

        """
        dnla20 = 1e-20 * nd_plasma_electron_line

        # ========================================================================

        # ITER-1996 H-mode power threshold database
        # Fit to 1996 H-mode power threshold database: nominal

        # i_l_h_threshold = 1
        iterdd = self.calculate_iter1996_nominal(
            dnla20, b_plasma_toroidal_on_axis, rmajor
        )

        # Fit to 1996 H-mode power threshold database: upper bound
        # i_l_h_threshold = 2
        iterdd_ub = self.calculate_iter1996_upper(
            dnla20, b_plasma_toroidal_on_axis, rmajor
        )

        # Fit to 1996 H-mode power threshold database: lower bound
        # i_l_h_threshold = 3
        iterdd_lb = self.calculate_iter1996_lower(
            dnla20, b_plasma_toroidal_on_axis, rmajor
        )

        # ========================================================================

        # Snipes 1997 ITER H-mode power threshold

        # i_l_h_threshold = 4
        snipes_1997 = self.calculate_snipes1997_iter(
            dnla20, b_plasma_toroidal_on_axis, rmajor
        )

        # i_l_h_threshold = 5
        snipes_1997_kappa = self.calculate_snipes1997_kappa(
            dnla20, b_plasma_toroidal_on_axis, rmajor, kappa
        )

        # ========================================================================

        # Martin et al (2008) for recent ITER scaling, with mass correction
        # and 95% confidence limits

        # i_l_h_threshold = 6
        martin_nominal = self.calculate_martin08_nominal(
            dnla20, b_plasma_toroidal_on_axis, a_plasma_surface, m_ions_total_amu
        )

        # i_l_h_threshold = 7
        martin_ub = self.calculate_martin08_upper(
            dnla20, b_plasma_toroidal_on_axis, a_plasma_surface, m_ions_total_amu
        )

        # i_l_h_threshold = 8
        martin_lb = self.calculate_martin08_lower(
            dnla20, b_plasma_toroidal_on_axis, a_plasma_surface, m_ions_total_amu
        )

        # ========================================================================

        # Snipes et al (2000) scaling with mass correction
        # Nominal, upper and lower

        # i_l_h_threshold = 9
        snipes_2000 = self.calculate_snipes2000_nominal(
            dnla20, b_plasma_toroidal_on_axis, rmajor, rminor, m_ions_total_amu
        )

        # i_l_h_threshold = 10
        snipes_2000_ub = self.calculate_snipes2000_upper(
            dnla20, b_plasma_toroidal_on_axis, rmajor, rminor, m_ions_total_amu
        )

        # i_l_h_threshold = 11
        snipes_2000_lb = self.calculate_snipes2000_lower(
            dnla20, b_plasma_toroidal_on_axis, rmajor, rminor, m_ions_total_amu
        )

        # ========================================================================

        # Snipes et al (2000) scaling (closed divertor) with mass correction
        # Nominal, upper and lower

        # i_l_h_threshold = 12
        snipes_2000_cd = self.calculate_snipes2000_closed_divertor_nominal(
            dnla20, b_plasma_toroidal_on_axis, rmajor, m_ions_total_amu
        )

        # i_l_h_threshold = 13
        snipes_2000_cd_ub = self.calculate_snipes2000_closed_divertor_upper(
            dnla20, b_plasma_toroidal_on_axis, rmajor, m_ions_total_amu
        )

        # i_l_h_threshold = 14
        snipes_2000_cd_lb = self.calculate_snipes2000_closed_divertor_lower(
            dnla20, b_plasma_toroidal_on_axis, rmajor, m_ions_total_amu
        )

        # ========================================================================

        # Hubbard et al. 2012 L-I threshold scaling

        # i_l_h_threshold = 15
        hubbard_2012 = self.calculate_hubbard2012_nominal(plasma_current, dnla20)

        # i_l_h_threshold = 16
        hubbard_2012_lb = self.calculate_hubbard2012_lower(plasma_current, dnla20)

        # i_l_h_threshold = 17
        hubbard_2012_ub = self.calculate_hubbard2012_upper(plasma_current, dnla20)

        # ========================================================================

        # Hubbard et al. 2017 L-I threshold scaling

        # i_l_h_threshold = 18
        hubbard_2017 = self.calculate_hubbard2017(
            dnla20, a_plasma_surface, b_plasma_toroidal_on_axis
        )

        # ========================================================================

        # Aspect ratio corrected Martin et al (2008)

        # i_l_h_threshold = 19
        martin_nominal_aspect = self.calculate_martin08_aspect_nominal(
            dnla20, b_plasma_toroidal_on_axis, a_plasma_surface, m_ions_total_amu, aspect
        )

        # i_l_h_threshold = 20
        martin_ub_aspect = self.calculate_martin08_aspect_upper(
            dnla20, b_plasma_toroidal_on_axis, a_plasma_surface, m_ions_total_amu, aspect
        )

        # i_l_h_threshold = 21
        martin_lb_aspect = self.calculate_martin08_aspect_lower(
            dnla20, b_plasma_toroidal_on_axis, a_plasma_surface, m_ions_total_amu, aspect
        )

        # ========================================================================

        return [
            iterdd,
            iterdd_ub,
            iterdd_lb,
            snipes_1997,
            snipes_1997_kappa,
            martin_nominal,
            martin_ub,
            martin_lb,
            snipes_2000,
            snipes_2000_ub,
            snipes_2000_lb,
            snipes_2000_cd,
            snipes_2000_cd_ub,
            snipes_2000_cd_lb,
            hubbard_2012,
            hubbard_2012_lb,
            hubbard_2012_ub,
            hubbard_2017,
            martin_nominal_aspect,
            martin_ub_aspect,
            martin_lb_aspect,
        ]

    def output_l_h_threshold_powers(self) -> None:
        """Output L-H transition power thresholds to file."""
        po.oheadr(self.outfile, "H-mode Power Threshold Scalings :")

        po.ovarin(
            self.outfile,
            "L-H threshold scaling switch",
            "(i_l_h_threshold)",
            physics_variables.i_l_h_threshold,
        )
        po.ocmmnt(
            self.outfile,
            f"{PlasmaConfinementTransitionModel(physics_variables.i_l_h_threshold).full_name}",
        )
        if (numerics.ioptimz > 0) and (numerics.active_constraints[14]):
            po.ovarre(
                self.outfile,
                "L-H threshold power (MW)",
                "(p_l_h_threshold_mw)",
                physics_variables.p_l_h_threshold_mw,
                "OP ",
            )
        else:
            po.ovarre(
                self.outfile,
                "L-H threshold power (NOT enforced) (MW)",
                "(p_l_h_threshold_mw)",
                physics_variables.p_l_h_threshold_mw,
                "OP ",
            )
        po.oblnkl(self.outfile)

        po.ovarre(
            self.outfile,
            "ITER 1996 scaling: nominal (MW)",
            "(l_h_threshold_powers(1))",
            physics_variables.l_h_threshold_powers[0],
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "ITER 1996 scaling: upper bound (MW)",
            "(l_h_threshold_powers(2))",
            physics_variables.l_h_threshold_powers[1],
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "ITER 1996 scaling: lower bound (MW)",
            "(l_h_threshold_powers(3))",
            physics_variables.l_h_threshold_powers[2],
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "ITER 1997 scaling (1) (MW)",
            "(l_h_threshold_powers(4))",
            physics_variables.l_h_threshold_powers[3],
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "ITER 1997 scaling (2) (MW)",
            "(l_h_threshold_powers(5))",
            physics_variables.l_h_threshold_powers[4],
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Martin 2008 scaling: nominal (MW)",
            "(l_h_threshold_powers(6))",
            physics_variables.l_h_threshold_powers[5],
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Martin 2008 scaling: 95% upper bound (MW)",
            "(l_h_threshold_powers(7))",
            physics_variables.l_h_threshold_powers[6],
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Martin 2008 scaling: 95% lower bound (MW)",
            "(l_h_threshold_powers(8))",
            physics_variables.l_h_threshold_powers[7],
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Snipes 2000 scaling: nominal (MW)",
            "(l_h_threshold_powers(9))",
            physics_variables.l_h_threshold_powers[8],
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Snipes 2000 scaling: upper bound (MW)",
            "(l_h_threshold_powers(10))",
            physics_variables.l_h_threshold_powers[9],
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Snipes 2000 scaling: lower bound (MW)",
            "(l_h_threshold_powers(11))",
            physics_variables.l_h_threshold_powers[10],
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Snipes 2000 scaling (closed divertor): nominal (MW)",
            "(l_h_threshold_powers(12))",
            physics_variables.l_h_threshold_powers[11],
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Snipes 2000 scaling (closed divertor): upper bound (MW)",
            "(l_h_threshold_powers(13))",
            physics_variables.l_h_threshold_powers[12],
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Snipes 2000 scaling (closed divertor): lower bound (MW)",
            "(l_h_threshold_powers(14))",
            physics_variables.l_h_threshold_powers[13],
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Hubbard 2012 L-I threshold - nominal (MW)",
            "(l_h_threshold_powers(15))",
            physics_variables.l_h_threshold_powers[14],
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Hubbard 2012 L-I threshold - lower bound (MW)",
            "(l_h_threshold_powers(16))",
            physics_variables.l_h_threshold_powers[15],
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Hubbard 2012 L-I threshold - upper bound (MW)",
            "(l_h_threshold_powers(17))",
            physics_variables.l_h_threshold_powers[16],
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Hubbard 2017 L-I threshold",
            "(l_h_threshold_powers(18))",
            physics_variables.l_h_threshold_powers[17],
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Martin 2008 aspect ratio corrected scaling: nominal (MW)",
            "(l_h_threshold_powers(19))",
            physics_variables.l_h_threshold_powers[18],
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Martin 2008 aspect ratio corrected scaling: 95% upper bound (MW)",
            "(l_h_threshold_powers(20))",
            physics_variables.l_h_threshold_powers[19],
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Martin 2008 aspect ratio corrected scaling: 95% lower bound (MW)",
            "(l_h_threshold_powers(21))",
            physics_variables.l_h_threshold_powers[20],
            "OP ",
        )
        po.oblnkl(self.outfile)
        if physics_variables.i_l_h_threshold in {9, 10, 11}:
            if (physics_variables.b_plasma_toroidal_on_axis < 0.78e0) or (
                physics_variables.b_plasma_toroidal_on_axis > 7.94e0
            ):
                po.ocmmnt(
                    self.outfile,
                    "(b_plasma_toroidal_on_axis outside Snipes 2000 fitted range)",
                )
                logger.warning(
                    "b_plasma_toroidal_on_axis outside Snipes 2000 fitted range"
                )

            if (physics_variables.rminor < 0.15e0) or (
                physics_variables.rminor > 1.15e0
            ):
                po.ocmmnt(self.outfile, "(rminor outside Snipes 2000 fitted range)")
                logger.warning("rminor outside Snipes 2000 fitted range")

            if (physics_variables.rmajor < 0.55e0) or (
                physics_variables.rmajor > 3.37e0
            ):
                po.ocmmnt(
                    self.outfile,
                    "(physics_variables.rmajor outside Snipes 2000 fitted range)",
                )
                logger.warning("rmajor outside Snipes 2000 fitted range")

            if (physics_variables.nd_plasma_electron_line < 0.09e20) or (
                physics_variables.nd_plasma_electron_line > 3.16e20
            ):
                po.ocmmnt(
                    self.outfile,
                    "(physics_variables.nd_plasma_electron_line outside Snipes 2000 "
                    "fitted range)",
                )
                logger.warning(
                    "nd_plasma_electron_line outside Snipes 2000 fitted range"
                )

            if (physics_variables.kappa < 1.0e0) or (physics_variables.kappa > 2.04e0):
                po.ocmmnt(
                    self.outfile,
                    "(physics_variables.kappa outside Snipes 2000 fitted range)",
                )
                logger.warning("kappa outside Snipes 2000 fitted range")

            if (physics_variables.triang < 0.07e0) or (
                physics_variables.triang > 0.74e0
            ):
                po.ocmmnt(self.outfile, "(triang outside Snipes 2000 fitted range)")
                logger.warning("triang outside Snipes 2000 fitted range")

        if physics_variables.i_l_h_threshold in {12, 13, 14}:
            po.ocmmnt(
                self.outfile,
                "(L-H threshold for closed divertor only. Limited data used in Snipes "
                "fit)",
            )
            po.oblnkl(self.outfile)
            logger.warning("Closed divertor only. Limited data used in Snipes fit")
        po.oblnkl(self.outfile)

    @staticmethod
    def calculate_iter1996_nominal(
        dnla20: float, b_plasma_toroidal_on_axis: float, rmajor: float
    ) -> float:
        """Calculate the nominal ITER-1996 L-H transition power threshold.

        Parameters
        ----------
        dnla20 : float
            Line averaged electron density in units of 10^20 m^-3.
        b_plasma_toroidal_on_axis : float
            Toroidal magnetic field [T]
        rmajor : float
            Plasma major radius [m]

        Returns
        -------
        float
            The ITER-1996 L-H transition power threshold [MW]

        References
        ----------
            - T. Takizuka and International Atomic Energy Agency, Vienna (Austria),
            "Threshold power and energy confinement for ITER". 1996.

            - J. C. Wesley, “International Thermonuclear Experimental Reactor: Physics
            issues, capabilities and physics program plans,”
            Physics of Plasmas, vol. 4, no. 7, pp. 2642-2652, Jul. 1997,
            doi: https://doi.org/10.1063/1.872406.
        """
        return 0.45 * dnla20**0.75 * b_plasma_toroidal_on_axis * rmajor**2

    @staticmethod
    def calculate_iter1996_upper(
        dnla20: float, b_plasma_toroidal_on_axis: float, rmajor: float
    ) -> float:
        """Calculate the upper variant ITER-1996 L-H transition power threshold.

        Parameters
        ----------
        dnla20 : float
            Line averaged electron density in units of 10^20 m^-3.
        b_plasma_toroidal_on_axis : float
            Toroidal magnetic field [T]
        rmajor : float
            Plasma major radius [m]

        Returns
        -------
        float
            The ITER-1996 L-H transition power threshold [MW]

        References
        ----------
            - T. Takizuka and International Atomic Energy Agency, Vienna (Austria),
            "Threshold power and energy confinement for ITER". 1996.

            - J. C. Wesley, “International Thermonuclear Experimental Reactor: Physics
            issues, capabilities and physics program plans,”
            Physics of Plasmas, vol. 4, no. 7, pp. 2642-2652, Jul. 1997,
            doi: https://doi.org/10.1063/1.872406.
        """
        return 0.3960502816 * dnla20 * b_plasma_toroidal_on_axis * rmajor**2.5

    @staticmethod
    def calculate_iter1996_lower(
        dnla20: float, b_plasma_toroidal_on_axis: float, rmajor: float
    ) -> float:
        """Calculate the lower variant ITER-1996 L-H transition power threshold.

        Parameters
        ----------
        dnla20 : float
            Line averaged electron density in units of 10^20 m^-3.
        b_plasma_toroidal_on_axis : float
            Toroidal magnetic field [T]
        rmajor : float
            Plasma major radius [m]

        Returns
        -------
        float
            The ITER-1996 L-H transition power threshold [MW]

        References
        ----------
            - T. Takizuka and International Atomic Energy Agency, Vienna (Austria),
            "Threshold power and energy confinement for ITER". 1996.

            - J. C. Wesley, “International Thermonuclear Experimental Reactor: Physics
            issues, capabilities and physics program plans,”
            Physics of Plasmas, vol. 4, no. 7, pp. 2642-2652, Jul. 1997,
            doi: https://doi.org/10.1063/1.872406.
        """
        return 0.5112987149 * dnla20**0.5 * b_plasma_toroidal_on_axis * rmajor**1.5

    @staticmethod
    def calculate_snipes1997_iter(
        dnla20: float, b_plasma_toroidal_on_axis: float, rmajor: float
    ) -> float:
        """Calculate the Snipes 1997 ITER L-H transition power threshold.

        Parameters
        ----------
        dnla20 : float
            Line averaged electron density in units of 10^20 m^-3.
        b_plasma_toroidal_on_axis : float
            Toroidal magnetic field [T]
        rmajor : float
            Plasma major radius [m]

        Returns
        -------
        float
            The Snipes 1997 L-H transition power threshold [MW]

        References
        ----------
            - J. A. Snipes and the ITER H-mode Threshold Database Working Group,
            "An Analysis of the H-mode Threshold in ITER,"
            Controlled Fusion and Plasma Physics, 24th EPS Conference,
            Berchtesgaden, June 9th-13th 1997, vol.21A, part III, p.961.
            url:https://library.ipp.mpg.de/EPS_24_Vol3_1997.pdf.
            *This is a conference poster*
        """
        return 0.65 * dnla20**0.93 * b_plasma_toroidal_on_axis**0.86 * rmajor**2.15

    @staticmethod
    def calculate_snipes1997_kappa(
        dnla20: float, b_plasma_toroidal_on_axis: float, rmajor: float, kappa: float
    ) -> float:
        """Calculate the Snipes 1997 ITER L-H transition power threshold with kappa
        factor.

        Parameters
        ----------
        dnla20 : float
            Line averaged electron density in units of 10^20 m^-3.
        b_plasma_toroidal_on_axis : float
            Toroidal magnetic field [T]
        rmajor : float
            Plasma major radius [m]
        kappa : float
            Plasma elongation

        Returns
        -------
        float
            The Snipes 1997 L-H transition power threshold with kappa factor [MW]

        References
        ----------
            - J. A. Snipes and the ITER H-mode Threshold Database Working Group,
            "An Analysis of the H-mode Threshold in ITER,"
            Controlled Fusion and Plasma Physics, 24th EPS Conference,
            Berchtesgaden, June 9th-13th 1997, vol.21A, part III, p.961.
            url:https://library.ipp.mpg.de/EPS_24_Vol3_1997.pdf.
            *This is a conference poster*
        """
        return (
            0.42
            * dnla20**0.80
            * b_plasma_toroidal_on_axis**0.90
            * rmajor**1.99
            * kappa**0.76
        )

    @staticmethod
    def calculate_martin08_nominal(
        dnla20: float,
        b_plasma_toroidal_on_axis: float,
        a_plasma_surface: float,
        m_ions_total_amu: float,
    ) -> float:
        """Calculate the nominal Martin L-H transition power threshold.

        Parameters
        ----------
        dnla20 : float
            Line averaged electron density in units of 10^20 m^-3.
        b_plasma_toroidal_on_axis : float
            Toroidal magnetic field [T]
        a_plasma_surface : float
            Plasma surface area [m^2]
        m_ions_total_amu : float
            Total ion mass in atomic mass units [amu]

        Returns
        -------
        float
            The Martin L-H transition power threshold [MW]

        Notes
        -----
            - A scaling with the total ion mass is used in this model. Martin 08 shows
            that P_LH scales with 1/m_i. It is stated; "When this mass dependence is
            applied to the deuterium-tritium discharges for ITER, the above predicted
            values of P_LH can be reduced by ~ 20%". We thus apply a (2/m_i) addition so
            that for a 50/50 D-T mixture (M_i = 2.5 amu), the predicted values is
            20% lower.


        References
        ----------
            - Y. R. Martin, T. Takizuka, and the I. C. H-mode. T. D. Group,
            “Power requirement for accessing the H-mode in ITER,”
            Journal of Physics: Conference Series, vol. 123, p. 012033, Jul. 2008,
            doi: https://doi.org/10.1088/1742-6596/123/1/012033.
        """
        return (
            0.0488
            * dnla20**0.717
            * b_plasma_toroidal_on_axis**0.803
            * a_plasma_surface**0.941
            * (2.0 / m_ions_total_amu)
        )

    @staticmethod
    def calculate_martin08_upper(
        dnla20: float,
        b_plasma_toroidal_on_axis: float,
        a_plasma_surface: float,
        m_ions_total_amu: float,
    ) -> float:
        """Calculate the upper Martin L-H transition power threshold.

        Parameters
        ----------
        dnla20 : float
            Line averaged electron density in units of 10^20 m^-3.
        b_plasma_toroidal_on_axis : float
            Toroidal magnetic field [T]
        a_plasma_surface : float
            Plasma surface area [m^2]
        m_ions_total_amu : float
            Total ion mass in atomic mass units [amu]

        Returns
        -------
        float
            The Martin L-H transition power threshold [MW]

        Notes
        -----
            - A scaling with the total ion mass is used in this model. Martin 08 shows
            that P_LH scales with 1/m_i. It is stated; "When this mass dependence is
            applied to the deuterium-tritium discharges for ITER, the above predicted
            values of P_LH can be reduced by ~ 20%". We thus apply a (2/m_i) addition
            so that for a 50/50 D-T mixture (M_i = 2.5 amu), the predicted values is
            20% lower.

        References
        ----------
            - Y. R. Martin, T. Takizuka, and the I. C. H-mode. T. D. Group, “Power
            requirement for accessing the H-mode in ITER,”
            Journal of Physics: Conference Series, vol. 123, p. 012033, Jul. 2008,
            doi: https://doi.org/10.1088/1742-6596/123/1/012033.
        """
        return (
            0.05166240355
            * dnla20**0.752
            * b_plasma_toroidal_on_axis**0.835
            * a_plasma_surface**0.96
            * (2.0 / m_ions_total_amu)
        )

    @staticmethod
    def calculate_martin08_lower(
        dnla20: float,
        b_plasma_toroidal_on_axis: float,
        a_plasma_surface: float,
        m_ions_total_amu: float,
    ) -> float:
        """Calculate the lower Martin L-H transition power threshold.

        Parameters
        ----------
        dnla20 : float
            Line averaged electron density in units of 10^20 m^-3.
        b_plasma_toroidal_on_axis : float
            Toroidal magnetic field [T]
        a_plasma_surface : float
            Plasma surface area [m^2]
        m_ions_total_amu : float
            Total ion mass in atomic mass units [amu]

        Returns
        -------
        float
            The Martin L-H transition power threshold [MW]

        Notes
        -----
            - A scaling with the total ion mass is used in this model. Martin 08 shows
            that P_LH scales with 1/m_i. It is stated; "When this mass dependence is
            applied to the deuterium-tritium discharges for ITER, the above predicted
            values of P_LH can be reduced by ~ 20%". We thus apply a (2/m_i) addition
            so that for a 50/50 D-T mixture (M_i = 2.5 amu), the predicted values is
            20% lower.

        References
        ----------
            - Y. R. Martin, T. Takizuka, and the I. C. H-mode. T. D. Group,
            “Power requirement for accessing the H-mode in ITER,”
            Journal of Physics: Conference Series, vol. 123, p. 012033, Jul. 2008,
            doi: https://doi.org/10.1088/1742-6596/123/1/012033.
        """
        return (
            0.04609619059
            * dnla20**0.682
            * b_plasma_toroidal_on_axis**0.771
            * a_plasma_surface**0.922
            * (2.0 / m_ions_total_amu)
        )

    @staticmethod
    def calculate_snipes2000_nominal(
        dnla20: float,
        b_plasma_toroidal_on_axis: float,
        rmajor: float,
        rminor: float,
        m_ions_total_amu: float,
    ) -> float:
        """Calculate the nominal Snipes 2000 L-H transition power threshold.

        Parameters
        ----------
        dnla20 : float
            Line averaged electron density in units of 10^20 m^-3.
        b_plasma_toroidal_on_axis : float
            Toroidal magnetic field [T]
        rmajor : float
            Plasma major radius [m]
        rminor : float
            Plasma minor radius [m]
        m_ions_total_amu : float
            Total ion mass in atomic mass units [amu]

        Returns
        -------
        float
            The Snipes 2000 L-H transition power threshold [MW]

        Notes
        -----
            - A scaling with the total ion mass is used in this model. Snipes cites
            that P_LH scales with 1/m_i. It is stated; "This results in a 20% reduction
            in the threshold power for a 50/50 D-T mixture compared with the pure
            deuterium results above". We thus apply a (2/m_i) addition so that for a
            50/50 D-T mixture (M_i = 2.5 amu), the predicted values is 20% lower.

        References
        ----------
            - J. A. Snipes and the I. H-mode. T. Group, “Latest results on the H-mode
            threshold using the international H-mode threshold database,”
            Plasma Physics and Controlled Fusion, vol. 42, no. 5A, pp. A299-A308,
            May 2000, doi: https://doi.org/10.1088/0741-3335/42/5a/336.
        """
        return (
            1.42
            * dnla20**0.58
            * b_plasma_toroidal_on_axis**0.82
            * rmajor
            * rminor**0.81
            * (2.0 / m_ions_total_amu)
        )

    @staticmethod
    def calculate_snipes2000_upper(
        dnla20: float,
        b_plasma_toroidal_on_axis: float,
        rmajor: float,
        rminor: float,
        m_ions_total_amu: float,
    ) -> float:
        """Calculate the upper Snipes 2000 L-H transition power threshold.

        Parameters
        ----------
        dnla20 : float
            Line averaged electron density in units of 10^20 m^-3.
        b_plasma_toroidal_on_axis : float
            Toroidal magnetic field [T]
        rmajor : float
            Plasma major radius [m]
        rminor : float
            Plasma minor radius [m]
        m_ions_total_amu : float
            Total ion mass in atomic mass units [amu]

        Returns
        -------
        float
            The Snipes 2000 L-H transition power threshold [MW]

        Notes
        -----
            - A scaling with the total ion mass is used in this model. Snipes cites
            that P_LH scales with 1/m_i. It is stated; "This results in a 20% reduction
            in the threshold power for a 50/50 D-T mixture compared with the pure
            deuterium results above". We thus apply a (2/m_i) addition so that for a
            50/50 D-T mixture (M_i = 2.5 amu), the predicted values is 20% lower.

        References
        ----------
            - J. A. Snipes and the I. H-mode. T. Group, “Latest results on the H-mode
            threshold using the international H-mode threshold database,”
            Plasma Physics and Controlled Fusion, vol. 42, no. 5A, pp. A299-A308,
            May 2000, doi: https://doi.org/10.1088/0741-3335/42/5a/336.

        """
        return (
            1.547
            * dnla20**0.615
            * b_plasma_toroidal_on_axis**0.851
            * rmajor**1.089
            * rminor**0.876
            * (2.0 / m_ions_total_amu)
        )

    @staticmethod
    def calculate_snipes2000_lower(
        dnla20: float,
        b_plasma_toroidal_on_axis: float,
        rmajor: float,
        rminor: float,
        m_ions_total_amu: float,
    ) -> float:
        """Calculate the lower Snipes 2000 L-H transition power threshold.

        Parameters
        ----------
        dnla20 : float
            Line averaged electron density in units of 10^20 m^-3.
        b_plasma_toroidal_on_axis : float
            Toroidal magnetic field [T]
        rmajor : float
            Plasma major radius [m]
        rminor : float
            Plasma minor radius [m]
        m_ions_total_amu : float
            Total ion mass in atomic mass units [amu]

        Returns
        -------
        float
            The Snipes 2000 L-H transition power threshold [MW]

        Notes
        -----
            - A scaling with the total ion mass is used in this model. Snipes cites
            that P_LH scales with 1/m_i. It is stated; "This results in a 20% reduction
            in the threshold power for a 50/50 D-T mixture compared with the pure
            deuterium results above". We thus apply a (2/m_i) addition so that for a
            50/50 D-T mixture (M_i = 2.5 amu), the predicted values is 20% lower.

        References
        ----------
            - J. A. Snipes and the I. H-mode. T. Group, “Latest results on the H-mode
            threshold using the international H-mode threshold database,”
            Plasma Physics and Controlled Fusion, vol. 42, no. 5A, pp. A299-A308,
            May 2000, doi: https://doi.org/10.1088/0741-3335/42/5a/336.

        """
        return (
            1.293
            * dnla20**0.545
            * b_plasma_toroidal_on_axis**0.789
            * rmajor**0.911
            * rminor**0.744
            * (2.0 / m_ions_total_amu)
        )

    @staticmethod
    def calculate_snipes2000_closed_divertor_nominal(
        dnla20: float,
        b_plasma_toroidal_on_axis: float,
        rmajor: float,
        m_ions_total_amu: float,
    ) -> float:
        """Calculate the nominal Snipes 2000 Closed Divertor L-H transition power
        threshold with CD factor.

        Parameters
        ----------
        dnla20 : float
            Line averaged electron density in units of 10^20 m^-3.
        b_plasma_toroidal_on_axis : float
            Toroidal magnetic field [T]
        rmajor : float
            Plasma major radius [m]
        m_ions_total_amu : float
            Total ion mass in atomic mass units [amu]

        Returns
        -------
        float
            The Snipes 2000 L-H transition power threshold with CD factor [MW]

        Notes
        -----
            - A scaling with the total ion mass is used in this model. Snipes cites
            that P_LH scales with 1/m_i. It is stated;m"This results in a 20% reduction
            in the threshold power for a 50/50 D-T mixture compared with the pure
            deuterium results above". We thus apply a (2/m_i) addition so that for a
            50/50 D-T mixture (M_i = 2.5 amu), the predicted values is 20% lower.

        References
        ----------
            - J. A. Snipes and the I. H-mode. T. Group, “Latest results on the H-mode
            threshold using the international H-mode threshold database,”
            Plasma Physics and Controlled Fusion, vol. 42, no. 5A, pp. A299-A308,
            May 2000, doi: https://doi.org/10.1088/0741-3335/42/5a/336.

        """
        return (
            0.8
            * dnla20**0.5
            * b_plasma_toroidal_on_axis**0.53
            * rmajor**1.51
            * (2.0 / m_ions_total_amu)
        )

    @staticmethod
    def calculate_snipes2000_closed_divertor_upper(
        dnla20: float,
        b_plasma_toroidal_on_axis: float,
        rmajor: float,
        m_ions_total_amu: float,
    ) -> float:
        """Calculate the upper Snipes 2000 Closed Divertor L-H transition power
        threshold with CD factor.

        Parameters
        ----------
        dnla20 : float
            Line averaged electron density in units of 10^20 m^-3.
        b_plasma_toroidal_on_axis : float
            Toroidal magnetic field [T]
        rmajor : float
            Plasma major radius [m]
        m_ions_total_amu : float
            Total ion mass in atomic mass units [amu]

        Returns
        -------
        float
            The Snipes 2000 L-H transition power threshold with CD factor [MW]

        Notes
        -----
            - A scaling with the total ion mass is used in this model. Snipes cites that
            P_LH scales with 1/m_i. It is stated; "This results in a 20% reduction in
            the threshold power for a 50/50 D-T mixture compared with the pure deuterium
            results above". We thus apply a (2/m_i) addition so that for a 50/50 D-T
            mixture (M_i = 2.5 amu), the predicted values is 20% lower.

        References
        ----------
            - J. A. Snipes and the I. H-mode. T. Group, “Latest results on the H-mode
            threshold using the international H-mode threshold database,”
            Plasma Physics and Controlled Fusion, vol. 42, no. 5A, pp. A299-A308,
            May 2000, doi: https://doi.org/10.1088/0741-3335/42/5a/336.

        """
        return (
            0.867
            * dnla20**0.561
            * b_plasma_toroidal_on_axis**0.588
            * rmajor**1.587
            * (2.0 / m_ions_total_amu)
        )

    @staticmethod
    def calculate_snipes2000_closed_divertor_lower(
        dnla20: float,
        b_plasma_toroidal_on_axis: float,
        rmajor: float,
        m_ions_total_amu: float,
    ) -> float:
        """Calculate the lower Snipes 2000 Closed Divertor L-H transition power
        threshold with CD factor.

        Parameters
        ----------
        dnla20 : float
            Line averaged electron density in units of 10^20 m^-3.
        b_plasma_toroidal_on_axis : float
            Toroidal magnetic field [T]
        rmajor : float
            Plasma major radius [m]
        m_ions_total_amu : float
            Total ion mass in atomic mass units [amu]

        Returns
        -------
        float
            The Snipes 2000 L-H transition power threshold with CD factor [MW]

        Notes
        -----
            - A scaling with the total ion mass is used in this model. Snipes cites that
              P_LH scales with 1/m_i. It is stated; "This results in a 20% reduction in
              the threshold power for a 50/50 D-T mixture compared with the pure
              deuterium results above". We thus apply a (2/m_i) addition so that for a
              50/50 D-T mixture (M_i = 2.5 amu), the predicted values is 20% lower.

        References
        ----------
            - J. A. Snipes and the I. H-mode. T. Group, “Latest results on the H-mode
            threshold using the international H-mode threshold database,”
            Plasma Physics and Controlled Fusion, vol. 42, no. 5A, pp. A299-A308,
            May 2000, doi: https://doi.org/10.1088/0741-3335/42/5a/336.

        """
        return (
            0.733
            * dnla20**0.439
            * b_plasma_toroidal_on_axis**0.472
            * rmajor**1.433
            * (2.0 / m_ions_total_amu)
        )

    @staticmethod
    def calculate_hubbard2012_nominal(plasma_current: float, dnla20: float) -> float:
        """Calculate the nominal Hubbard 2012 L-I transition power threshold.

        Parameters
        ----------
        plasma_current : float
            Plasma current [A]
        dnla20 : float
            Line averaged electron density in units of 10^20 m^-3.

        Returns
        -------
        float
            The Hubbard 2012 L-I transition power threshold [MW]

        References
        ----------
        - A. E. Hubbard et al., “Threshold conditions for transitions to I-mode and
          H-mode with unfavourable ion grad B drift direction,”
          Nuclear Fusion, vol. 52, no. 11, pp. 114009-114009, Oct. 2012,
          doi: https://doi.org/10.1088/0029-5515/52/11/114009.

        """
        return 2.11 * (plasma_current / 1e6) ** 0.94 * dnla20**0.65

    @staticmethod
    def calculate_hubbard2012_upper(plasma_current: float, dnla20: float) -> float:
        """Calculate the upper Hubbard 2012 L-I transition power threshold.

        Parameters
        ----------
        plasma_current : float
            Plasma current [A]
        dnla20 : float
            Line averaged electron density in units of 10^20 m^-3.

        Returns
        -------
        float
            The Hubbard 2012 L-I transition power threshold [MW]

        References
        ----------
        - A. E. Hubbard et al., “Threshold conditions for transitions to I-mode and
          H-mode with unfavourable ion grad B drift direction,”
          Nuclear Fusion, vol. 52, no. 11, pp. 114009-114009, Oct. 2012,
          doi: https://doi.org/10.1088/0029-5515/52/11/114009.

        """
        return 2.11 * (plasma_current / 1e6) ** 1.18 * dnla20**0.83

    @staticmethod
    def calculate_hubbard2012_lower(plasma_current: float, dnla20: float) -> float:
        """Calculate the lower Hubbard 2012 L-I transition power threshold.

        Parameters
        ----------
        plasma_current : float
            Plasma current [A]
        dnla20 : float
            Line averaged electron density in units of 10^20 m^-3.

        Returns
        -------
        float
            The Hubbard 2012 L-I transition power threshold [MW]

        References
        ----------
        - A. E. Hubbard et al., “Threshold conditions for transitions to I-mode and
          H-mode with unfavourable ion grad B drift direction,”
          Nuclear Fusion, vol. 52, no. 11, pp. 114009-114009, Oct. 2012,
          doi: https://doi.org/10.1088/0029-5515/52/11/114009.

        """
        return 2.11 * (plasma_current / 1e6) ** 0.7 * dnla20**0.47

    @staticmethod
    def calculate_hubbard2017(
        dnla20: float, a_plasma_surface: float, b_plasma_toroidal_on_axis: float
    ) -> float:
        """Calculate the Hubbard 2017 L-I transition power threshold.

        Parameters
        ----------
        dnla20 : float
            Line averaged electron density in units of 10^20 m^-3.
        a_plasma_surface : float
            Plasma surface area [m^2]
        b_plasma_toroidal_on_axis : float
            Toroidal magnetic field [T]

        Returns
        -------
        float
            The Hubbard 2017 L-I transition power threshold [MW]

        Notes
        -----
        - The scaling is given in the caption of Figure 6 in the reference.

        References
        ----------
        - A. E. Hubbard et al., “Physics and performance of the I-mode regime over an
          expanded operating space on Alcator C-Mod,” Nuclear Fusion, vol. 57, no. 12,
          p. 126039, Oct. 2017, doi: https://doi.org/10.1088/1741-4326/aa8570.

        """
        return 0.162 * dnla20 * a_plasma_surface * b_plasma_toroidal_on_axis**0.26

    @staticmethod
    def calculate_martin08_aspect_nominal(
        dnla20: float,
        b_plasma_toroidal_on_axis: float,
        a_plasma_surface: float,
        m_ions_total_amu: float,
        aspect: float,
    ) -> float:
        """Calculate the nominal Martin L-H transition power threshold with aspect ratio
        correction from T Takizuka.

        Parameters
        ----------
        dnla20 : float
            Line averaged electron density in units of 10^20 m^-3.
        b_plasma_toroidal_on_axis : float
            Toroidal magnetic field [T]
        a_plasma_surface : float
            Plasma surface area [m^2]
        m_ions_total_amu : float
            Total ion mass in atomic mass units [amu]
        aspect : float
            Plasma aspect ratio

        Returns
        -------
        float
            The Martin L-H transition power threshold [MW]

        Notes
        -----
            - Thus will return an aspect ratio correction of the aspect ratio is less
            than or equal to 2.7. If not the usual Martin 2008 scaling will be returned.

            - A scaling with the total ion mass is used in this model. Martin 08 shows
            that P_LH scales with 1/m_i. It is stated; "When this mass dependence is
            applied to the deuterium-tritium discharges for ITER, the above predicted
            values of P_LH can be reduced by ~ 20%". We thus apply a (2/m_i) addition
            so that for a 50/50 D-T mixture (M_i = 2.5 amu), the predicted values is
            20% lower.

        References
        ----------
            - Y. R. Martin, T. Takizuka, and the I. C. H-mode. T. D. Group, “Power
            requirement for accessing the H-mode in ITER,”
            Journal of Physics: Conference Series, vol. 123, p. 012033, Jul. 2008,
            doi: https://doi.org/10.1088/1742-6596/123/1/012033.

            - T. Takizuka et.al, “Roles of aspect ratio, absolute B and effective Z of
            the H-mode power threshold in tokamaks of the ITPA database,”
            Plasma Physics and Controlled Fusion, vol. 46, no. 5A, pp. A227-A233,
            Apr. 2004, doi: https://doi.org/10.1088/0741-3335/46/5a/024.

        """
        if aspect <= 2.7:
            aspect_correction = 0.098 * aspect / (1.0 - (2.0 / (1.0 + aspect)) ** 0.5)
        else:
            aspect_correction = 1.0

        return (
            0.0488
            * dnla20**0.717
            * b_plasma_toroidal_on_axis**0.803
            * a_plasma_surface**0.941
            * (2.0 / m_ions_total_amu)
            * aspect_correction
        )

    @staticmethod
    def calculate_martin08_aspect_upper(
        dnla20: float,
        b_plasma_toroidal_on_axis: float,
        a_plasma_surface: float,
        m_ions_total_amu: float,
        aspect: float,
    ) -> float:
        """Calculate the upper Martin L-H transition power threshold with aspect ratio
        correction from T Takizuka.

        Parameters
        ----------
        dnla20 : float
            Line averaged electron density in units of 10^20 m^-3.
        b_plasma_toroidal_on_axis : float
            Toroidal magnetic field [T]
        a_plasma_surface : float
            Plasma surface area [m^2]
        m_ions_total_amu : float
            Total ion mass in atomic mass units [amu]
        aspect : float
            Plasma aspect ratio

        Returns
        -------
        float
            The Martin L-H transition power threshold [MW]

        Notes
        -----
            - Thus will return an aspect ratio correction of the aspect ratio is less
            than or equal to 2.7. If not the usual Martin 2008 scaling will be returned.

            - A scaling with the total ion mass is used in this model. Martin 08 shows
            that P_LH scales with 1/m_i. It is stated; "When this mass dependence is
            applied to the deuterium-tritium discharges for ITER, the above predicted
            values of P_LH can be reduced by ~ 20%". We thus apply a (2/m_i) addition
            so that for a 50/50 D-T mixture (M_i = 2.5 amu), the predicted values is
            20% lower.

        References
        ----------
            - Y. R. Martin, T. Takizuka, and the I. C. H-mode. T. D. Group, “Power
            requirement for accessing the H-mode in ITER,”
            Journal of Physics: Conference Series, vol. 123, p. 012033, Jul. 2008,
            doi: https://doi.org/10.1088/1742-6596/123/1/012033.

            - T. Takizuka et.al, “Roles of aspect ratio, absolute B and effective Z of
            the H-mode power threshold in tokamaks of the ITPA database,”
            Plasma Physics and Controlled Fusion, vol. 46, no. 5A, pp. A227-A233,
            Apr. 2004, doi: https://doi.org/10.1088/0741-3335/46/5a/024.

        """
        if aspect <= 2.7:
            aspect_correction = 0.098 * aspect / (1.0 - (2.0 / (1.0 + aspect)) ** 0.5)
        else:
            aspect_correction = 1.0

        return (
            0.05166240355
            * dnla20**0.752
            * b_plasma_toroidal_on_axis**0.835
            * a_plasma_surface**0.96
            * (2.0 / m_ions_total_amu)
            * aspect_correction
        )

    @staticmethod
    def calculate_martin08_aspect_lower(
        dnla20: float,
        b_plasma_toroidal_on_axis: float,
        a_plasma_surface: float,
        m_ions_total_amu: float,
        aspect: float,
    ) -> float:
        """Calculate the lower Martin L-H transition power threshold with aspect ratio
        correction from T Takizuka.

        Parameters
        ----------
        dnla20 : float
            Line averaged electron density in units of 10^20 m^-3.
        b_plasma_toroidal_on_axis : float
            Toroidal magnetic field [T]
        a_plasma_surface : float
            Plasma surface area [m^2]
        m_ions_total_amu : float
            Total ion mass in atomic mass units [amu]
        aspect : float
            Plasma aspect ratio

        Returns
        -------
        float
            The Martin L-H transition power threshold [MW]

        Notes
        -----
            - Thus will return an aspect ratio correction of the aspect ratio is less
            than or equal to 2.7. if not the usual Martin 2008 scaling will be returned.

            - A scaling with the total ion mass is used in this model. Martin 08 shows
            that P_LH scales with 1/m_i. It is stated; "When this mass dependence is
            applied to the deuterium-tritium discharges for ITER, the above predicted
            values of P_LH can be reduced by ~ 20%". We thus apply a (2/m_i) addition
            so that for a 50/50 D-T mixture (M_i = 2.5 amu), the predicted values is
            20% lower.

        References
        ----------
            - Y. R. Martin, T. Takizuka, and the I. C. H-mode. T. D. Group, “Power
            requirement for accessing the H-mode in ITER,”
            Journal of Physics: Conference Series, vol. 123, p. 012033, Jul. 2008,
            doi: https://doi.org/10.1088/1742-6596/123/1/012033.

            - T. Takizuka et.al, “Roles of aspect ratio, absolute B and effective Z of
            the H-mode power threshold in tokamaks of the ITPA database,”
            Plasma Physics and Controlled Fusion, vol. 46, no. 5A, pp. A227-A233,
            Apr. 2004, doi: https://doi.org/10.1088/0741-3335/46/5a/024.

        """
        if aspect <= 2.7:
            aspect_correction = 0.098 * aspect / (1.0 - (2.0 / (1.0 + aspect)) ** 0.5)
        else:
            aspect_correction = 1.0

        return (
            0.04609619059
            * dnla20**0.682
            * b_plasma_toroidal_on_axis**0.771
            * a_plasma_surface**0.922
            * (2.0 / m_ions_total_amu)
            * aspect_correction
        )
