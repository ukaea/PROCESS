import logging
from dataclasses import dataclass

import matplotlib.pyplot as plt
import numpy as np
import numpy.typing as npt
import scipy.constants as const
from scipy.interpolate import interp1d
from skimage import measure

import process.io.mfile as mf

logger = logging.getLogger(__name__)


@dataclass
class ExtremalPoint:
    elongation: float
    triangularity: float
    x_point: bool
    squareness: float = 0.0

    def __post_init__(self):
        self.elongation = float(self.elongation)
        self.triangularity = float(self.triangularity)
        self.x_point = bool(self.x_point)
        self.squareness = float(self.squareness)

        if self.elongation <= 0:
            raise ValueError(f"Elongation must be positive: {self.elongation}")
        if abs(self.triangularity) > 1:
            raise ValueError(
                f"Triangularity must be between -1 and 1: {self.triangularity}"
            )
        if abs(self.squareness) > 0.5:
            raise ValueError(
                f"squareness must be between -0.5 and 0.5: {self.squareness}"
            )

    @classmethod
    def at_coordinates(
        cls,
        r: float,
        z: float,
        x_point: bool,
        eps: float,
        rmajor: float,
    ):
        elongation = abs(z) / rmajor / eps
        triangularity = (rmajor - r) / rmajor / eps

        if x_point:
            elongation /= 1.1
            triangularity /= 1.1

        return cls(elongation, triangularity, x_point)

    @property
    def get_x_point(self) -> bool:
        return self.x_point


class AnalyticGradShafranovSolution:
    """

        Class to represent an analytic Grad-Shafranov solution for an up-down axisymmetric plasma equilibrium.

        :references:
        A. J. Cerfon and J. P. Freidberg, “One size fits all analytic solutions to the Grad-Shafranov equation,”
        vol. 17, no. 3, pp. 032502-032502, Mar. 2010,
        doi: https://doi.org/10.1063/1.3328818.
    ‌
    """

    __slots__ = (
        "b_plasma_toroidal_on_axis",
        "beta_normalised",
        "beta_poloidal",
        "beta_toroidal",
        "beta_total",
        "boundary_height",
        "boundary_radius",
        "c_plasma_anticlockwise",
        "c_plasma_ma",
        "coefficients",
        "eps",
        "equatorial_point_inner_xy",
        "equatorial_point_outer_xy",
        "kink_safety_factor",
        "lower_point",
        "lower_point_xy",
        "magnetic_axis",
        "normalised_circumference",
        "normalised_volume",
        "poloidal_to_toroidal_flux",
        "pressure_parameter",
        "psi_0",
        "psi_axis",
        "psi_separatrix_toroidal",
        "q_profile",
        "rmajor",
        "toroidal_field_anticlockwise",
        "toroidal_to_poloidal_flux",
        "upper_point",
        "upper_point_xy",
    )

    def __init__(
        self,
        rmajor: float,
        pressure_parameter: float,
        eps: float,
        upper_point: ExtremalPoint,
        lower_point: ExtremalPoint,
        b_plasma_toroidal_on_axis: float,
        c_plasma_ma: float,
        kink_safety_factor: float | None = None,
        c_plasma_anticlockwise: bool = True,
        toroidal_field_anticlockwise: bool = True,
        use_d_shaped_model: bool = False,
    ):
        """
        Parameters
        ----------
        rmajor: float
            Major radius of plasma [m].
        pressure_parameter: float
            Parameter controlling the size of the pressure function A.
            Larger A gives a higher pressure. This can only really be set
            via trial and error to get the desired beta.
        eps: float
            Inverse aspect ratio epsilon [] = minor radius [m] / major radius [m].
        elongation: float
            Plasma elongation kappa []. Can be a tuple with the upper and lower triangularity.
        triangularity: float
            Plasma triangularity delta []. Can be a tuple with the upper and lower elongation.
        b_plasma_toroidal_on_axis: float
            Magnetic field strength at the geometric axis r = rmajor [T].
        c_plasma_ma: float
            Total plasma current [MA].
        kink_safety_factor: float, optional
            Kink safety factor q_star. If None (default), this is calculated using the plasma current.
            Otherwise the value of the plasma current is calculated using the provided value.
        """
        self.rmajor: float = float(rmajor)
        self.pressure_parameter: float = float(pressure_parameter)
        self.eps: float = float(eps)

        if not isinstance(upper_point, ExtremalPoint):
            raise ValueError("upper_point must be ExtremalPoint")
        if not isinstance(lower_point, ExtremalPoint):
            raise ValueError("lower_point must be ExtremalPoint")

        self.upper_point = upper_point
        self.lower_point = lower_point

        self.b_plasma_toroidal_on_axis: float = abs(float(b_plasma_toroidal_on_axis))
        self.c_plasma_anticlockwise: bool = c_plasma_anticlockwise
        self.toroidal_field_anticlockwise: bool = toroidal_field_anticlockwise

        if not self.toroidal_field_anticlockwise:
            self.b_plasma_toroidal_on_axis *= -1

        # Solve for the weighting coefficients for each of the psi polynomials.
        self.calculate_polynomial_psi_coefficients()

        # Calculate magnetic axis location.
        self.calculate_magnetic_axis()

        # Calculate (r, z) of boundary contour.
        self.calcuate_boundary_contour()

        # Calculate the normalised circumference and volume.
        self.calculate_geometry_factors(use_d_shaped_model=use_d_shaped_model)

        # Use either the plasma current or kink safety factor to calculate the other.
        e, b_plasma_toroidal_on_axis = self.eps, self.b_plasma_toroidal_on_axis
        rmajor, len_plasma_poloidal_norm = self.rmajor, self.normalised_circumference

        if kink_safety_factor is None:
            self.c_plasma_ma: float = abs(float(c_plasma_ma))
            self.kink_safety_factor = abs(
                e
                * b_plasma_toroidal_on_axis
                * rmajor
                * len_plasma_poloidal_norm
                / const.mu_0
                / (1e6 * self.c_plasma_ma)
            )
        else:
            self.kink_safety_factor = abs(float(kink_safety_factor))
            self.c_plasma_ma = 1e-6 * abs(
                e
                * b_plasma_toroidal_on_axis
                * rmajor
                * len_plasma_poloidal_norm
                / const.mu_0
                / self.kink_safety_factor
            )

        if not self.c_plasma_anticlockwise:
            self.c_plasma_ma *= -1

        # Set dummy value of psi axis. This will be set in calculate_metrics() to match the prescribed plasma current.
        self.psi_0 = 1.0
        self.calculate_metrics()
        self.calculate_q_profile()

        # Add interpolator to convert poloidal flux to toroidal flux.
        self.add_poloidal_toroidal_convertor()

    @property
    def upper_elongation(self) -> float:
        return self.upper_point.elongation

    @property
    def lower_elongation(self) -> float:
        return self.lower_point.elongation

    @property
    def upper_triangularity(self) -> float:
        return self.upper_point.triangularity

    @property
    def lower_triangularity(self) -> float:
        return self.lower_point.triangularity

    @property
    def upper_squareness(self) -> float:
        return self.upper_point.squareness

    @property
    def lower_squareness(self) -> float:
        return self.lower_point.squareness

    @staticmethod
    def psi_polynomials(x: float, y: float) -> tuple[float]:
        """
        The set of homogeneous polynomials (psi terms) that solve the Grad-Shafranov equation.

        :Notes: Polynomial numbering is given in Equation 8 of Cerfon and Freidberg (2010).
        """
        psi_1 = 1
        psi_2 = x**2
        psi_3 = y**2 - x**2 * np.log(x)
        psi_4 = x**2 * (x**2 - 4 * y**2)
        psi_5 = y**2 * (2 * y**2 - 9 * x**2) + x**2 * np.log(x) * (3 * x**2 - 12 * y**2)
        psi_6 = x**2 * (x**4 - 12 * x**2 * y**2 + 8 * y**4)
        psi_7 = (-15 * x**4 + 180 * y**2 * x**2 - 120 * y**4) * x**2 * np.log(x) + (
            75 * x**4 - 140 * x**2 * y**2 + 8 * y**4
        ) * y**2
        psi_8 = y
        psi_9 = y * x**2
        psi_10 = y * (y**2 - 3 * x**2 * np.log(x))
        psi_11 = y * (3 * x**4 - 4 * y**2 * x**2)
        psi_12 = y * (
            (8 * y**2 - 80 * x**2 * np.log(x)) * y**2 + (60 * np.log(x) - 45) * x**4
        )

        return (
            psi_1,
            psi_2,
            psi_3,
            psi_4,
            psi_5,
            psi_6,
            psi_7,
            psi_8,
            psi_9,
            psi_10,
            psi_11,
            psi_12,
        )

    @staticmethod
    def psi_polynomials_dx(x: float, y: float) -> tuple[float]:
        """
        First derivative with respect to x of homogeneous polynomials (psi terms) that
        solve the Grad-Shafranov equation.

        :Notes: Polynomial derivatives not given in reference. Manually derived.
        """

        dpsi_1_dx = 0
        dpsi_2_dx = 2 * x
        dpsi_3_dx = -x * (1 + 2 * np.log(x))
        dpsi_4_dx = x * (4 * x**2 - 8 * y**2)
        dpsi_5_dx = 3 * x * ((4 * x**2 - 8 * y**2) * np.log(x) + x**2 - 10 * y**2)
        dpsi_6_dx = x * (6 * x**4 - 48 * x**2 * y**2 + 16 * y**4)
        dpsi_7_dx = (
            -5
            * x
            * (
                (18 * x**4 - 144 * x**2 * y**2 + 48 * y**4) * np.log(x)
                + 3 * x**4
                - 96 * x**2 * y**2
                + 80 * y**4
            )
        )
        dpsi_8_dx = 0
        dpsi_9_dx = 2 * x * y
        dpsi_10_dx = -3 * x * y * (2 * np.log(x) + 1)
        dpsi_11_dx = x * y * (12 * x**2 - 8 * y**2)
        dpsi_12_dx = (
            40 * x * y * ((6 * x**2 - 4 * y**2) * np.log(x) - 3 * x**2 - 2 * y**2)
        )
        return (
            dpsi_1_dx,
            dpsi_2_dx,
            dpsi_3_dx,
            dpsi_4_dx,
            dpsi_5_dx,
            dpsi_6_dx,
            dpsi_7_dx,
            dpsi_8_dx,
            dpsi_9_dx,
            dpsi_10_dx,
            dpsi_11_dx,
            dpsi_12_dx,
        )

    @staticmethod
    def psi_polynomials_dy(x: float, y: float) -> tuple[float]:
        """
        First derivative with respect to y of homogeneous polynomials (psi terms) that
        solve the Grad-Shafranov equation.

        :Notes: Polynomial derivatives not given in reference. Manually derived.
        """
        dpsi_1_dy = 0
        dpsi_2_dy = 0
        dpsi_3_dy = 2 * y
        dpsi_4_dy = -8 * (x**2) * y
        dpsi_5_dy = y * (8 * (y**2) - (x**2) * (18 + 24 * np.log(x)))
        dpsi_6_dy = y * (-24 * (x**4) + 32 * (x**2) * (y**2))
        dpsi_7_dy = y * (
            48 * (y**4)
            - (480 * np.log(x) + 560) * (x**2) * (y**2)
            + (360 * np.log(x) + 150) * (x**4)
        )
        dpsi_8_dy = 1
        dpsi_9_dy = x**2
        dpsi_10_dy = 3 * ((y**2) - (x**2) * np.log(x))
        dpsi_11_dy = (x**2) * (3 * (x**2) - 12 * (y**2))
        dpsi_12_dy = 40 * (y**4) + 15 * (x**2) * (
            ((-16 * (y**2) + 4 * (x**2)) * np.log(x)) - 3 * (x**2)
        )
        return (
            dpsi_1_dy,
            dpsi_2_dy,
            dpsi_3_dy,
            dpsi_4_dy,
            dpsi_5_dy,
            dpsi_6_dy,
            dpsi_7_dy,
            dpsi_8_dy,
            dpsi_9_dy,
            dpsi_10_dy,
            dpsi_11_dy,
            dpsi_12_dy,
        )

    @staticmethod
    def psi_polynomials_dx2(x: float, y: float) -> tuple[float]:
        """
        Second derivative with respect to x of homogeneous polynomials (psi terms) that
        solve the Grad-Shafranov equation.

        :Notes: Polynomial derivatives not given in reference. Manually derived.
        """

        d2_psi_1_dx2 = 0
        d2_psi_2_dx2 = 2
        d2_psi_3_dx2 = -3 - 2 * np.log(x)
        d2_psi_4_dx2 = 12 * x**2 - 8 * y**2
        d2_psi_5_dx2 = (36 * x**2 - 24 * y**2) * np.log(x) + 21 * x**2 - 54 * y**2
        d2_psi_6_dx2 = 30 * x**4 - 144 * x**2 * y**2 + 16 * y**4
        d2_psi_7_dx2 = (
            (-450 * x**4 + 2160 * x**2 * y**2 - 240 * y**4) * np.log(x)
            - 165 * x**4
            + 2160 * x**2 * y**2
            - 640 * y**4
        )
        d2_psi_8_dx2 = 0
        d2_psi_9_dx2 = 2 * y
        d2_psi_10_dx2 = -3 * y * (2 * np.log(x) + 3)
        d2_psi_11_dx2 = y * (36 * x**2 - 8 * y**2)
        d2_psi_12_dx2 = y * (
            (720 * x**2 - 160 * y**2) * np.log(x) - 120 * x**2 - 240 * y**2
        )
        return (
            d2_psi_1_dx2,
            d2_psi_2_dx2,
            d2_psi_3_dx2,
            d2_psi_4_dx2,
            d2_psi_5_dx2,
            d2_psi_6_dx2,
            d2_psi_7_dx2,
            d2_psi_8_dx2,
            d2_psi_9_dx2,
            d2_psi_10_dx2,
            d2_psi_11_dx2,
            d2_psi_12_dx2,
        )

    @staticmethod
    def psi_polynomials_dy2(x: float, y: float) -> tuple[float]:
        """
        Second derivative with respect to y of homogeneous polynomials (psi terms) that
        solve the Grad-Shafranov equation.

        :Notes: Polynomial derivatives not given in reference. Manually derived.
        """
        d2_psi_1_dy2 = 0
        d2_psi_2_dy2 = 0
        d2_psi_3_dy2 = 2
        d2_psi_4_dy2 = -8 * x**2
        d2_psi_5_dy2 = 24 * y**2 - x**2 * (18 + 24 * np.log(x))
        d2_psi_6_dy2 = x**2 * (-24 * x**2 + 96 * y**2)
        d2_psi_7_dy2 = (
            240 * y**4
            - (1440 * np.log(x) + 1680) * x**2 * y**2
            + (360 * np.log(x) + 150) * x**4
        )
        d2_psi_8_dy2 = 0
        d2_psi_9_dy2 = 0
        d2_psi_10_dy2 = 6 * y
        d2_psi_11_dy2 = -24 * x**2 * y
        d2_psi_12_dy2 = y * (160 * y**2 - 480 * x**2 * np.log(x))
        return (
            d2_psi_1_dy2,
            d2_psi_2_dy2,
            d2_psi_3_dy2,
            d2_psi_4_dy2,
            d2_psi_5_dy2,
            d2_psi_6_dy2,
            d2_psi_7_dy2,
            d2_psi_8_dy2,
            d2_psi_9_dy2,
            d2_psi_10_dy2,
            d2_psi_11_dy2,
            d2_psi_12_dy2,
        )

    def psi_particular(self, x: float, y: float) -> float:
        """
        Particular solution of the normalised Grad-Shafranov equation. x = r / rmajor and y = z / rmajor are
        the radius r and height z normalised to the major radius rmajor.

        :Notes: Equation 6 of Cerfon and Freidberg (2010).
        """
        return 0.5 * self.pressure_parameter * x**2 * np.log(x) - (x**4 / 8.0) * (
            1.0 + self.pressure_parameter
        )

    def psi_particular_dx(self, x: float, y: float) -> float:
        """First derivative of the particular solution with respect to x.

        :Notes: Manually derived.
        """
        return (
            0.5
            * x
            * (
                self.pressure_parameter * (2 * np.log(x) + 1)
                - x**2 * (1 + self.pressure_parameter)
            )
        )

    def psi_particular_dx2(self, x: float, y: float) -> float:
        """Second derivative of the particular solution with respect to x.
        :Notes: Manually derived.
        """

        return (
            self.pressure_parameter * (1.5 + np.log(x))
            - 1.5 * (1.0 + self.pressure_parameter) * x**2
        )

    def psi_bar(self, x: float, y: float) -> float:
        """
        Poloidal flux function normalised to psi0. This is NOT the commonly encountered psi normalised psiN!
        psi_bar is zero at the separatrix and some positive value at the magnetic axis.
        """
        psi = self.psi_particular(x, y)

        # Add weighted sum of the homogeneous solutions.
        for c, p in zip(self.coefficients, self.psi_polynomials(x, y), strict=False):
            psi += c * p

        return psi

    def psi_bar_dx(self, x: float, y: float) -> float:
        """First derivative of poloidal flux function normalised to psi0 with respect to x."""
        psi_dx = self.psi_particular_dx(x, y)

        # Add weighted sum of the homogeneous solutions.
        for c, p in zip(self.coefficients, self.psi_polynomials_dx(x, y), strict=False):
            psi_dx += c * p

        return psi_dx

    def psi_bar_dy(self, x: float, y: float) -> float:
        """First derivative of poloidal flux function normalised to psi0 with respect to y."""
        psi_dy = 0

        # Add weighted sum of the homogeneous solutions.
        for c, p in zip(self.coefficients, self.psi_polynomials_dy(x, y), strict=False):
            psi_dy += c * p

        return psi_dy

    def psi_bar_dx2(self, x: float, y: float) -> float:
        """Second derivative of poloidal flux function normalised to psi0 with respect to x."""
        psi_dx2 = self.psi_particular_dx2(x, y)

        # Add weighted sum of the homogeneous solutions.
        for c, p in zip(self.coefficients, self.psi_polynomials_dx2(x, y), strict=False):
            psi_dx2 += c * p

        return psi_dx2

    def psi_bar_dy2(self, x: float, y: float) -> float:
        """Second derivative of poloidal flux function normalised to psi0 with respect to y."""
        psi_dy2 = 0

        # Add weighted sum of the homogeneous solutions.
        for c, p in zip(self.coefficients, self.psi_polynomials_dy2(x, y), strict=False):
            psi_dy2 += c * p

        return psi_dy2

    def psi(self, r: float, z: float) -> float:
        """Poloidal flux function [Wb]."""
        rmajor = self.rmajor
        return self.psi_0 * self.psi_bar(r / rmajor, z / rmajor)

    def psi_dr(self, r: float, z: float) -> float:
        """First derivative of the poloidal flux function with respect to r [Wb / m]."""
        rmajor = self.rmajor
        return self.psi_0 * self.psi_bar_dx(r / rmajor, z / rmajor) / rmajor

    def psi_dz(self, r: float, z: float) -> float:
        """First derivative of the poloidal flux function with respect to z [Wb / m]."""
        rmajor = self.rmajor
        return self.psi_0 * self.psi_bar_dy(r / rmajor, z / rmajor) / rmajor

    def psi_norm(self, r: float, z: float) -> float:
        rmajor = self.rmajor
        return self.psi_bar_to_psi_norm(self.psi_bar(r / rmajor, z / rmajor))

    def magnetic_field(self, r: float, z: float) -> tuple[float, float, float]:
        """(r, phi, z) components of the magnetic field [T]."""
        psi_norm = self.psi_norm(r, z)
        # NOTE: Need a test for this!
        b_r = self.psi_dz(r, z) / r
        b_z = -self.psi_dr(r, z) / r

        # We account for the direction of the toroidal field in f_function.
        b_toroidal = self.f_function(psi_norm) / r
        # If plasma current is negative poloidal field reverses sign.
        if not self.c_plasma_anticlockwise:
            b_r *= -1

        return b_r, b_toroidal, b_z

    def calculate_polynomial_psi_coefficients(self):
        """
        Solve for the weighted coefficients of the polynomials defining Ψ. We fit to a d shaped contour with the
        required geometry factors at 3 points:
            Inner equatorial point: point of minimum r at midplane (z=0) on the boundary contour.
            Outer equatorial point: point of minimum r at midplane (z=0) on the boundary contour.
            High point: point of maximum z on the boundary contour.
            Upper X point: point of maximum z on the boundary contour.

        :Notes: Equations given as Equations 10, 11 and 12 in Cerfon and Freidberg (2010).

        """
        e = self.eps

        # Some coefficients from D shaped model. Use average elongation, triangularity and squareness at midplane.

        k_mid = 0.5 * (self.upper_elongation + self.lower_elongation)
        d_mid = 0.5 * (self.upper_triangularity + self.lower_triangularity)
        s_mid = 0.5 * (self.upper_squareness + self.lower_squareness)
        alpha_mid = np.arcsin(d_mid)
        n1_mid = -((1 + alpha_mid) ** 2) / self.eps / k_mid**2 / (1 + 2 * s_mid**2)
        n2_mid = (1 - alpha_mid) ** 2 / self.eps / k_mid**2 / (1 + 2 * s_mid**2)

        # Points to fit D shaped model at.
        self.equatorial_point_inner_xy = (1 - self.eps, 0)
        self.equatorial_point_outer_xy = (1 + self.eps, 0)

        # We solve the system y = Mx to find the coefficient vector x.
        matrix = np.zeros((12, 12))
        y = np.zeros(12)

        # Outer equatorial point (Ψ = 0).
        matrix[0] = self.psi_polynomials(*self.equatorial_point_outer_xy)
        y[0] = -self.psi_particular(*self.equatorial_point_outer_xy)

        # Inner equatorial point (Ψ = 0).
        matrix[1] = self.psi_polynomials(*self.equatorial_point_inner_xy)
        y[1] = -self.psi_particular(*self.equatorial_point_inner_xy)

        # Outer equatorial point up down symmetry (d(Ψ)/dy = 0).
        matrix[2] = self.psi_polynomials_dy(*self.equatorial_point_outer_xy)

        # Inner equatorial point up down symmetry (d(Ψ)/dy = 0).
        matrix[3] = self.psi_polynomials_dy(*self.equatorial_point_inner_xy)

        # Outer equatorial point curvature (d^2(Ψ)/dy^2 + N1 * d(Ψ)/dx = 0).
        matrix[4] = np.array(
            self.psi_polynomials_dy2(*self.equatorial_point_outer_xy)
        ) + n1_mid * np.array(self.psi_polynomials_dx(*self.equatorial_point_outer_xy))
        y[4] = -n1_mid * self.psi_particular_dx(*self.equatorial_point_outer_xy)

        # Inner equatorial point curvature (d^2(Ψ)/dy^2 + N2 * d(Ψ)/dx = 0).
        matrix[5] = np.array(
            self.psi_polynomials_dy2(*self.equatorial_point_inner_xy)
        ) + n2_mid * np.array(self.psi_polynomials_dx(*self.equatorial_point_inner_xy))
        y[5] = -n2_mid * self.psi_particular_dx(*self.equatorial_point_inner_xy)

        if self.upper_point.get_x_point:
            k, d = self.upper_elongation, self.upper_triangularity
            self.upper_point_xy = (1 - 1.1 * d * self.eps, 1.1 * k * self.eps)

            # Upper X point (Ψ = 0).
            matrix[6] = self.psi_polynomials(*self.upper_point_xy)
            y[6] = -self.psi_particular(*self.upper_point_xy)

            # B poloidal = 0 at upper X point (d(Ψ)/dx = 0).
            matrix[7] = self.psi_polynomials_dx(*self.upper_point_xy)
            y[7] = -self.psi_particular_dx(*self.upper_point_xy)

            # B poloidal = 0 at upper X point (d(Ψ)/dy = 0).
            matrix[8] = self.psi_polynomials_dy(*self.upper_point_xy)
        else:
            # Upper high point.
            k, d, s = (
                self.upper_elongation,
                self.upper_triangularity,
                self.upper_squareness,
            )
            n_3 = -k * (1 - 2 * s**2) / e / (1 - d**2)
            self.upper_point_xy = (1 - d * self.eps, k * self.eps)

            # Upper high point (Ψ = 0).
            matrix[6] = self.psi_polynomials(*self.upper_point_xy)
            y[6] = -self.psi_particular(*self.upper_point_xy)

            # Upper high point maximum (d(Ψ)/dx = 0).
            matrix[7] = self.psi_polynomials_dx(*self.upper_point_xy)
            y[7] = -self.psi_particular_dx(*self.upper_point_xy)

            # Upper high point curvature (d^2(Ψ)/dx^2 + n_3 * d(Ψ)/dy = 0).
            matrix[8] = np.array(
                self.psi_polynomials_dx2(*self.upper_point_xy)
            ) + n_3 * np.array(self.psi_polynomials_dy(*self.upper_point_xy))
            y[8] = -self.psi_particular_dx2(*self.upper_point_xy)

        if self.lower_point.get_x_point:
            k, d = self.lower_elongation, self.lower_triangularity
            self.lower_point_xy = (1 - 1.1 * d * self.eps, -1.1 * k * self.eps)

            # Lower X point (Ψ = 0).
            matrix[9] = self.psi_polynomials(*self.lower_point_xy)
            y[9] = -self.psi_particular(*self.lower_point_xy)

            # B poloidal = 0 at lower X point (d(Ψ)/dx = 0).
            matrix[10] = self.psi_polynomials_dx(*self.lower_point_xy)
            y[10] = -self.psi_particular_dx(*self.lower_point_xy)

            # B poloidal = 0 at lower X point (d(Ψ)/dy = 0).
            matrix[11] = self.psi_polynomials_dy(*self.lower_point_xy)
        else:
            # Lower high point.
            k, d, s = (
                self.lower_elongation,
                self.lower_triangularity,
                self.lower_squareness,
            )
            n_3 = k * (1 - 2 * s**2) / self.eps / (1 - d**2)
            self.lower_point_xy = (1 - d * self.eps, -k * self.eps)

            # Lower high point (Ψ = 0).
            matrix[9] = self.psi_polynomials(*self.lower_point_xy)
            y[9] = -self.psi_particular(*self.lower_point_xy)

            # Lower high point maximum (d(Ψ)/dx = 0).
            matrix[10] = self.psi_polynomials_dx(*self.lower_point_xy)
            y[10] = -self.psi_particular_dx(*self.lower_point_xy)

            # Lower high point curvature (d^2(Ψ)/dx^2 + n_3 * d(Ψ)/dy = 0).
            matrix[11] = np.array(
                self.psi_polynomials_dx2(*self.lower_point_xy)
            ) + n_3 * np.array(self.psi_polynomials_dy(*self.lower_point_xy))
            y[11] = -self.psi_particular_dx2(*self.lower_point_xy)

        self.coefficients = np.linalg.solve(matrix, y)

    def calculate_magnetic_axis(
        self, tolerance: float = 1.0e-6, max_iterations: int = 100
    ):
        """
        Calculate the magnetic axis using Newton's method.

        This method finds the location of the magnetic axis by solving for the point
        where the gradient of the poloidal flux (d(psi)/dx and d(psi)/dy) is zero.

        Args:
            tolerance (float): Convergence tolerance for the Newton's method.
            max_iterations (int): Maximum number of iterations allowed.

        Returns:
            None: Updates the `self.magnetic_axis` attribute with the calculated position.
        """
        # Initial guess for the magnetic axis (near the geometric center).
        magnetic_axis = np.array([1.0, 0.0])

        for i in range(max_iterations):
            # Compute the gradient and Hessian components of psi_bar at the current point.
            psi_dx = self.psi_bar_dx(*magnetic_axis)
            psi_dy = self.psi_bar_dy(*magnetic_axis)
            psi_dx2 = self.psi_bar_dx2(*magnetic_axis)
            psi_dy2 = self.psi_bar_dy2(*magnetic_axis)
            psi_dxy = (
                self.psi_bar_dx(magnetic_axis[0], magnetic_axis[1] + 1e-6) - psi_dx
            ) / 1e-6

            # Construct the Hessian matrix and gradient vector.
            hessian = np.array([[psi_dx2, psi_dxy], [psi_dxy, psi_dy2]])
            gradient = np.array([psi_dx, psi_dy])

            # Solve for the correction using the inverse of the Hessian.
            try:
                correction = np.linalg.solve(hessian, gradient)
            except np.linalg.LinAlgError:
                logger.error("Hessian is singular; cannot find magnetic axis.")
                break

            # Update the magnetic axis position.
            magnetic_axis -= correction

            # Check for convergence.
            if np.linalg.norm(correction) < tolerance:
                logger.info(f"Magnetic axis found in {i + 1} iterations.")
                self.magnetic_axis = magnetic_axis * self.rmajor
                return

        logger.error(
            "Failed to find magnetic axis within the maximum number of iterations."
        )

    @property
    def dr_shafranov(self) -> float:
        """Shafranov shift: radial distance between magnetic axis and geometric centre."""
        return self.magnetic_axis[0] - self.rmajor

    def d_shape_boundary(self, theta: float) -> float:
        """D shaped boundary contour for the prescribed geometry factors."""
        theta = np.array(theta)
        x, y = np.zeros_like(theta), np.zeros_like(theta)
        mask = np.logical_and(theta >= 0, theta <= np.pi)

        def d_shape(theta, k, alpha, s):
            x = 1 + self.eps * np.cos(theta + alpha * np.sin(theta))
            y = self.eps * k * np.sin(theta + s * np.sin(2 * theta))
            return x, y

        # Above midplane.
        k, d, s = self.upper_elongation, self.upper_triangularity, self.upper_squareness
        alpha = np.arcsin(d)
        x[mask], y[mask] = d_shape(theta[mask], k, alpha, s)

        # Below midplane.
        k, d, s = self.lower_elongation, self.lower_triangularity, self.lower_squareness
        alpha = np.arcsin(d)
        x[~mask], y[~mask] = d_shape(theta[~mask], k, alpha, s)

        return x, y

    def d_shape_boundary_derivatives(self, theta: float) -> float:
        theta = np.array(theta)
        xprime, yprime = np.zeros_like(theta), np.zeros_like(theta)
        mask = np.logical_and(theta >= 0, theta <= np.pi)

        e = self.eps

        def d_shape_prime(theta, k, alpha, s):
            xprime = (
                -e * np.sin(theta + alpha * np.sin(theta)) * (1 + alpha * np.cos(theta))
            )
            yprime = (
                e
                * k
                * np.cos(theta + s * np.sin(2 * theta))
                * (1 + 2 * s * np.cos(2 * theta))
            )
            return xprime, yprime

        # Above midplane.
        k, d, s = self.upper_elongation, self.upper_triangularity, self.upper_squareness
        alpha = np.arcsin(d)
        xprime[mask], yprime[mask] = d_shape_prime(theta[mask], k, alpha, s)

        # Below midplane.
        k, d, s = self.lower_elongation, self.lower_triangularity, self.lower_squareness
        alpha = np.arcsin(d)
        xprime[~mask], yprime[~mask] = d_shape_prime(theta[~mask], k, alpha, s)

        return xprime, yprime

    def calcuate_boundary_contour(
        self,
        n_points: int = 101,
        psi_norm_threshold: float = 1.0e-3,
    ):
        """Calculate boundary contour of psi map using Newton's method"""
        # Calculate the extremal value of psi so we can calculate psi norm.
        psi_bar_0 = self.psi_bar(*self.magnetic_axis / self.rmajor)

        # Distort d shaped value such that psi=0.
        theta = np.linspace(0, 2 * np.pi, n_points)
        x_d_shape, y_d_shape = self.d_shape_boundary(theta)
        x_boundary, y_boundary = np.copy(x_d_shape), np.copy(y_d_shape)

        for i, t in enumerate(theta):
            psi_bar = self.psi_bar(x_boundary[i], y_boundary[i])
            psi_norm = psi_bar / psi_bar_0

            # If psi norm is close enough to zero, break.
            if abs(psi_norm) < psi_norm_threshold:
                break

            # Use Netwon's method to update boundary position.
            dpsi_dx = self.psi_bar_dx(x_boundary[i], y_boundary[i])
            dpsi_dy = self.psi_bar_dy(x_boundary[i], y_boundary[i])

            cos_t, sin_t = np.cos(t), np.sin(t)
            dpsi_dv = cos_t * dpsi_dx + sin_t * dpsi_dy

            x_boundary[i] -= cos_t * psi_bar / dpsi_dv
            y_boundary[i] -= sin_t * psi_bar / dpsi_dv

        self.boundary_radius = x_boundary * self.rmajor
        self.boundary_height = y_boundary * self.rmajor

    def calculate_geometry_factors(self, use_d_shaped_model: bool = True):
        """
        Calculate the normalised circumference and volume of the plasma. Can either calculate based on the estimated
        boundary contour from the D shaped model or from the fitted poloidal flux function.
        """
        if use_d_shaped_model:
            theta = np.linspace(0, 2 * np.pi, 101)
            x, y = self.d_shape_boundary(theta)
            xprime, yprime = self.d_shape_boundary_derivatives(theta)
            rprime = (xprime**2 + yprime**2) ** 0.5

            self.normalised_circumference = np.trapz(rprime, theta)
            self.normalised_volume = -np.trapz(x * xprime * y, theta)
        else:
            circumference = 0
            volume = 0

            for dr, dz in zip(
                np.diff(self.boundary_radius),
                np.diff(self.boundary_height),
                strict=False,
            ):
                circumference += (dr**2 + dz**2) ** 0.5

            for i in range(len(self.boundary_radius) - 1):
                r_1, z_1 = self.boundary_radius[i], self.boundary_height[i]
                r_2, z_2 = self.boundary_radius[i + 1], self.boundary_height[i + 1]
                volume += 0.5 * (r_1 * z_1 + r_2 * z_2) * (r_2 - r_1)

            self.normalised_circumference = circumference / self.rmajor
            self.normalised_volume = -volume / self.rmajor**3

    def metric_computation_grid(self) -> tuple[npt.NDArray[float], npt.NDArray[float]]:
        """
        Grid of normalised radius x = r / rmajor and height y = z / rmajor used to calculate numerical integrals
        for calculating the plasma current and plasma beta.
        """
        e, k_up, k_down = (
            self.eps,
            self.upper_elongation,
            self.lower_elongation,
        )
        x = np.linspace(1 - e, 1 + e, 100)
        y = np.linspace(-e * k_down, e * k_up, 101)
        return x, y

    def calculate_metrics(self):
        """
        Calculate the poloidal flux coordinate normalisation (psi_0) and the plasma 'figures of merit':
        beta_poloidal, beta_toroidal, beta_total and beta_normalised.
        """

        len_plasma_poloidal_norm, vol_plasma_norm = (
            self.normalised_circumference,
            self.normalised_volume,
        )
        rmajor, b_plasma_toroidal_on_axis = self.rmajor, self.b_plasma_toroidal_on_axis
        pressure_parameter, qstar = self.pressure_parameter, self.kink_safety_factor

        # Calculate two integrals appearing in expressions for beta and plasma current.
        x, y = self.metric_computation_grid()
        x_grid, ygrid = np.meshgrid(x, y, indexing="ij")
        dxdy = (x[1] - x[0]) * (y[1] - y[0])

        # Ignore contribution from anything outside the separatrix.
        psi_bar = self.psi_bar(x_grid, ygrid)
        mask = psi_bar < 0

        # This is proportional to the total plasma current.
        ip_integrand = (pressure_parameter / x_grid) - (1 + pressure_parameter) * x_grid
        ip_integrand[mask] = 0
        ip_integral = np.sum(ip_integrand) * dxdy

        # This is proportional to the average plasma pressure.
        psix_integrand = psi_bar * x_grid
        psix_integrand[mask] = 0
        psix_integral = np.sum(psix_integrand) * dxdy

        # Calculate beta of plasma.
        self.beta_poloidal = (
            (
                2
                * len_plasma_poloidal_norm**2
                * (1 + pressure_parameter)
                / vol_plasma_norm
            )
            * psix_integral
            / ip_integral**2
        )
        self.beta_toroidal = self.beta_poloidal * self.eps**2 / qstar**2
        self.beta_total = self.beta_poloidal * self.eps**2 / (self.eps**2 + qstar**2)
        self.beta_normalised = (
            self.eps * rmajor * abs(b_plasma_toroidal_on_axis) / abs(self.c_plasma_ma)
        ) * self.beta_total

        # Calculate value of psi at magnetic axis for given plasma current.
        # Enforce psi to be increasing with minor radius.
        self.psi_0 = -abs(
            self.c_plasma_ma * 1e6 * const.mu_0 * self.rmajor / ip_integral
        )

        # Evaluate value of psi_norm at the magnetic axis.
        self.psi_axis = self.psi(*self.magnetic_axis)

    def calculate_q_profile(self, mesh_size: int = 30):
        """Calculate q profile."""
        rmajor = self.rmajor

        # The q profile q(psi) = 1/(2*pi) * int{F(psi) / (r |grad{psi}|)} dl_p where
        # l_p is the poloidal distance along the surface of constant psi.
        q_profile = np.zeros(mesh_size)

        # Calculate (x, y) locations of contours.
        xmesh, ymesh = self.metric_computation_grid()
        dxmesh, dymesh = xmesh[-1] - xmesh[0], ymesh[-1] - ymesh[0]
        psi_bar_axis = self.psi_bar(*self.magnetic_axis / rmajor)
        psi_bar_norm_grid = (
            self.psi_bar(*np.meshgrid(xmesh, ymesh, indexing="ij")) / psi_bar_axis
        )

        def integrand(x, y):
            # Return 1 / (r |grad{psi}|). Technically we use use
            # 1 / (x |grad_x,y{psi}|) but the factor of major radius cancels.
            dpsi_dx = self.psi_0 * self.psi_bar_dx(x, y)
            dpsi_dy = self.psi_0 * self.psi_bar_dy(x, y)
            mod_grad_psi = (dpsi_dx**2 + dpsi_dy**2) ** 0.5
            return 1 / (mod_grad_psi * x)

        def calculate_arclength(x, y):
            # Calculate arclength around contour.
            lp = np.zeros_like(x)
            for k, (dx, dy) in enumerate(zip(np.diff(x), np.diff(y), strict=False)):
                lp[k + 1] = lp[k] + (dx**2 + dy**2) ** 0.5

            return lp

        # NOTE: This is psi_bar normalised, so it is 1 at the magnetic axis
        # and zero at the boundary. We will use the computed boundary contour
        # to do psi_bar_norm = 0. Also skip psi_bar_norm = 1 as there is just a single point (bad numerics).
        psi_norm_mesh = np.linspace(1, 0, mesh_size + 1)[1:-1]

        for i, psi_norm in enumerate(psi_norm_mesh):
            # Use this instead of matplotlib.contour as the latter forces figure creation.
            contour = measure.find_contours(psi_bar_norm_grid, level=psi_norm)

            if len(contour) == 0:
                raise ValueError(f"Unable to find contour for psi norm = {psi_norm}")

            x = xmesh[0] + dxmesh * (contour[0][:, 0] / (psi_bar_norm_grid.shape[0] - 1))
            y = ymesh[0] + dymesh * (contour[0][:, 1] / (psi_bar_norm_grid.shape[1] - 1))

            # F function is a flux function so we can move it out the integral (F / r is toroidal field).
            # As we are COCOS 11 q > 0 so take absolute value of F.
            f = abs(self.f_function(psi_norm))

            # Integrate using trapezium rule.
            lp = rmajor * calculate_arclength(x, y)  # Poloidal arclength [m].
            q_profile[i] = f * np.trapz(integrand(x, y), lp)
        # Use pre-computed boundary contour for separatrix. As there is a
        # saddle point the matplotlib contours will sometimes follow the contours
        # towards the divertor instead of following the high field side boundary.
        x_bdy, y_bdy = self.boundary_radius / rmajor, self.boundary_height / rmajor

        # Integrate using trapezium rule. As we are COCOS 11 q > 0 so take absolute value of F.
        f = abs(self.f_function(1))
        lp = rmajor * calculate_arclength(x_bdy, y_bdy)  # Poloidal arclength [m].
        q_profile[-1] = f * np.trapz(integrand(x_bdy, y_bdy), lp)

        # Scale q profile by 2*pi to match definition of poloidal flux function psi.
        q_profile /= 2 * np.pi

        # Set q profile.
        self.q_profile = q_profile

    def add_poloidal_toroidal_convertor(self):
        """Define function to convert from poloidal flux to toroidal flux."""
        poloidal_flux = np.linspace(self.psi_axis, 0, len(self.q_profile))

        # Toroidal flux is int{q dpsi_poloidal}
        toroidal_flux = np.zeros_like(poloidal_flux)

        for i in range(len(toroidal_flux) - 1):
            dpsi_tor = (
                0.5
                * (self.q_profile[i] + self.q_profile[i + 1])
                * (poloidal_flux[i + 1] - poloidal_flux[i])
            )
            toroidal_flux[i + 1] = toroidal_flux[i] + dpsi_tor

        self.psi_separatrix_toroidal = toroidal_flux[-1]

        # Create some interpolators to convert between them.
        self.poloidal_to_toroidal_flux = interp1d(
            poloidal_flux,
            toroidal_flux,
            bounds_error=False,
            fill_value=(toroidal_flux[0], toroidal_flux[-1]),
        )
        self.toroidal_to_poloidal_flux = interp1d(
            toroidal_flux,
            poloidal_flux,
            bounds_error=False,
            fill_value=(poloidal_flux[0], poloidal_flux[-1]),
        )

    def psi_bar_to_psi_norm(self, psi_bar: float) -> float:
        """Convert psi_bar parameter used in the normalised Grad Shafranov equation to the standard normalised poloidal flux co-ordinate."""
        return 1 - (psi_bar * self.psi_0 / self.psi_axis)

    def psi_norm_to_psi_bar(self, psi_norm: float) -> float:
        """Convert normalised poloidal flux co-ordinate to the psi_bar parameter used in the normalised Grad Shafranov equation."""
        return (1 - psi_norm) * self.psi_axis / self.psi_0

    def psi_toroidal(self, r, z):
        """Toroidal flux function."""
        psi_poloidal = self.psi(r, z)
        return self.poloidal_to_toroidal_flux(psi_poloidal)

    def psi_norm_toroidal(self, r, z):
        psi_toroidal = self.psi_toroidal(r, z)
        return psi_toroidal / self.psi_separatrix_toroidal

    def psi_norm_poloidal_to_toroidal(self, psi_norm_poloidal):
        """Convert normalised poloidal flux coordinate to normalised toroidal flux coordinate."""
        psi_poloidal = self.psi_axis * (1 - psi_norm_poloidal)
        psi_toroidal = self.poloidal_to_toroidal_flux(psi_poloidal)
        return psi_toroidal / self.psi_separatrix_toroidal

    def psi_norm_toroidal_to_poloidal(self, psi_norm_toroidal):
        """Convert normalised toroidal flux coordinate to normalised poloidal flux coordinate."""
        psi_toroidal = self.psi_separatrix_toroidal * psi_norm_toroidal
        psi_poloidal = self.toroidal_to_poloidal_flux(psi_toroidal)
        return 1 - (psi_poloidal / self.psi_axis)

    def pressure_kpa(self, psi_norm: float):
        """Plasma pressure as a function of the normalised poloidal flux [kPa]."""
        pressure_parameter = self.pressure_parameter
        # Clip psi to 0 so there is not negative pressure.
        psi_bar = np.clip(self.psi_norm_to_psi_bar(psi_norm), 0, None)
        return (
            1e-3
            * (self.psi_0**2 / self.rmajor**4 / const.mu_0)
            * (1 + pressure_parameter)
            * psi_bar
        )

    def f_function(self, psi_norm: float):
        """F function radius * toroidal magnetic field as a function of the normalised poloidal flux [Wb/m]."""
        rmajor, b_plasma_toroidal_on_axis = self.rmajor, self.b_plasma_toroidal_on_axis
        psi0, pressure_parameter = self.psi_0, self.pressure_parameter
        # Clip psi to 0 to avoid unphysical magnetic fields.
        psi_bar = np.clip(self.psi_norm_to_psi_bar(psi_norm), 0, None)
        f = (
            rmajor
            * (
                b_plasma_toroidal_on_axis**2
                - (2 * psi0**2 / rmajor**4) * pressure_parameter * psi_bar
            )
            ** 0.5
        )

        # If toroidal field is clockwise flip the sign of f.
        if not self.toroidal_field_anticlockwise:
            f *= -1

        return f

    def toroidal_current_density_ka_per_m2(self, r) -> float:
        """Toroidal current density [kA m^-2]."""
        return (
            1e-3
            * self.psi_0
            * (
                self.pressure_parameter / (r / self.rmajor)
                - (1.0 + self.pressure_parameter) * (r / self.rmajor) ** 2
            )
            / (const.mu_0 * self.rmajor**3)
        )

    def toroidal_field(self, r: float, psi_norm: float) -> float:
        psi_bar = np.clip(self.psi_norm_to_psi_bar(psi_norm), 0, None)

        return np.sqrt(
            (self.rmajor**2 / (r / self.rmajor) ** 2)
            * (
                4.5**2
                - (2 * self.psi_0**2 / self.rmajor**4)
                * self.pressure_parameter
                * psi_bar
            )
        )

    def plotting_xy_arrays(self, padding=1.05, n_points: int = 100):
        """Grid of (x, y) points that encloses the entire plasma boundary plus some padding"""
        e = self.eps
        xmin, xmax = (1 - padding * e), (1 + padding * e)
        ymin, ymax = padding * self.lower_point_xy[1], padding * self.upper_point_xy[1]

        return np.linspace(xmin, xmax, n_points), np.linspace(ymin, ymax, n_points)

    def plotting_rz_arrays(self, **kwargs):
        """Arrays of (r, z) points that encloses the entire plasma boundary plus some padding."""
        x, y = self.plotting_xy_arrays(**kwargs)
        return x * self.rmajor, y * self.rmajor


class Limiter(AnalyticGradShafranovSolution):
    __slots__ = ()

    def __init__(
        self,
        rmajor: float,
        pressure_parameter: float,
        eps: float,
        elongation: float,
        triangularity: float,
        b_plasma_toroidal_on_axis: float,
        c_plasma_ma: float,
        kink_safety_factor,
        squareness: float = 0.0,
        c_plasma_anticlockwise: bool = True,
        toroidal_field_anticlockwise: bool = True,
        use_d_shaped_model: bool = False,
    ):
        upper_point = ExtremalPoint(
            elongation, triangularity, False, squareness=squareness
        )
        lower_point = upper_point
        super().__init__(
            rmajor,
            pressure_parameter,
            eps,
            upper_point,
            lower_point,
            b_plasma_toroidal_on_axis,
            c_plasma_ma,
            kink_safety_factor,
            c_plasma_anticlockwise,
            toroidal_field_anticlockwise,
            use_d_shaped_model,
        )


class DoubleNull(AnalyticGradShafranovSolution):
    __slots__ = ()

    def __init__(
        self,
        rmajor: float,
        pressure_parameter: float,
        eps: float,
        elongation: float,
        triangularity: float,
        b_plasma_toroidal_on_axis: float,
        c_plasma_ma: float,
        kink_safety_factor: float | None = None,
        squareness: float = 0.0,
        c_plasma_anticlockwise: bool = True,
        toroidal_field_anticlockwise: bool = True,
        use_d_shaped_model: bool = False,
    ):
        upper_point = ExtremalPoint(
            elongation, triangularity, True, squareness=squareness
        )
        lower_point = upper_point

        super().__init__(
            rmajor,
            pressure_parameter,
            eps,
            upper_point,
            lower_point,
            b_plasma_toroidal_on_axis,
            c_plasma_ma,
            kink_safety_factor,
            c_plasma_anticlockwise,
            toroidal_field_anticlockwise,
            use_d_shaped_model,
        )


class SingleNull(AnalyticGradShafranovSolution):
    __slots__ = ()

    def __init__(
        self,
        rmajor: float,
        pressure_parameter: float,
        eps: float,
        elongation: float,
        triangularity: float,
        b_plasma_toroidal_on_axis: float,
        c_plasma_ma: float,
        kink_safety_factor: float | None = None,
        squareness: float = 0.0,
        c_plasma_anticlockwise: bool = True,
        toroidal_field_anticlockwise: bool = True,
        use_d_shaped_model: bool = False,
        lower_x: bool = True,
    ):
        if lower_x:
            upper_point = ExtremalPoint(
                elongation, triangularity, False, squareness=squareness
            )
            lower_point = ExtremalPoint(
                elongation, triangularity, True, squareness=squareness
            )
        else:
            upper_point = ExtremalPoint(
                elongation, triangularity, True, squareness=squareness
            )
            lower_point = ExtremalPoint(
                elongation, triangularity, False, squareness=squareness
            )

        super().__init__(
            rmajor,
            pressure_parameter,
            eps,
            upper_point,
            lower_point,
            b_plasma_toroidal_on_axis,
            c_plasma_ma,
            kink_safety_factor,
            c_plasma_anticlockwise,
            toroidal_field_anticlockwise,
            use_d_shaped_model,
        )


def plot_analytic_equilibrium(
    axis: plt.Axes, mfile: mf.MFile, scan: int, fig: plt.Figure
):
    i_single_null = mfile.get("i_single_null", scan=scan)
    beta_thermal_toroidal_vol_avg = mfile.get("beta_thermal_toroidal_vol_avg", scan=scan)
    beta_thermal_poloidal_vol_avg = mfile.get("beta_thermal_poloidal_vol_avg", scan=scan)
    beta_thermal_vol_avg = mfile.get("beta_thermal_vol_avg", scan=scan)
    beta_norm_thermal = mfile.get("beta_norm_thermal", scan=scan)

    n_plasma_profile_elements = int(mfile.get("n_plasma_profile_elements", scan=scan))

    pres_plasma_thermal_total_profile = [
        mfile.get(f"pres_plasma_thermal_total_profile{i}", scan=scan)
        for i in range(n_plasma_profile_elements)
    ]
    rminor = mfile.get("rminor", scan=scan)

    # Helper to build plasma object for a given pressure_parameter
    def make_plasma(p_param):
        common_kwargs = {
            "rmajor": mfile.get("rmajor", scan=scan),
            "pressure_parameter": p_param,
            "elongation": mfile.get("kappa", scan=scan),
            "triangularity": mfile.get("triang", scan=scan),
            "b_plasma_toroidal_on_axis": mfile.get(
                "b_plasma_toroidal_on_axis", scan=scan
            ),
            "c_plasma_ma": mfile.get("plasma_current_ma", scan=scan),
            "kink_safety_factor": None,
        }
        if i_single_null == 1:
            # SingleNull expects eps expressed as 1/aspect in original code path
            common_kwargs["eps"] = 1.0 / mfile.get("aspect", scan=scan)
            # provide elongation fields explicitly if SingleNull needs separate upper/lower;
            # constructor used previously only passed 'elongation' and 'triangularity'
            return SingleNull(**common_kwargs)
        # DoubleNull uses eps directly from data in original code
        common_kwargs["eps"] = 1.0 / mfile.get("aspect", scan=scan)
        return DoubleNull(**common_kwargs)

    # Solve for pressure_parameter such that plasma.beta_toroidal == beta_thermal_toroidal_vol_avg
    target = float(beta_thermal_toroidal_vol_avg)

    # If target is NaN or zero-ish, skip solving and use a modest default
    if not np.isfinite(target) or abs(target) < 1e-15:
        p_opt = 0.2  # fallback default used previously
    else:

        def beta_diff(p):
            try:
                pl = make_plasma(p)
                return float(pl.beta_toroidal) - target
            except (ValueError, TypeError, ArithmeticError):
                # If construction fails for this p, return large positive diff to force bracket expansion
                return np.inf

        # Bracket search
        low = 1e-8
        high = 1.0
        f_low = beta_diff(low)
        f_high = beta_diff(high)

        # Expand high until we bracket a sign change or reach a max
        max_expand = 50
        expand_iter = 0
        while (
            not (np.isfinite(f_low) and np.isfinite(f_high) and f_low * f_high <= 0)
            and expand_iter < max_expand
        ):
            if not np.isfinite(f_low):
                low = max(low * 0.1, 1e-12)
                f_low = beta_diff(low)
            if not np.isfinite(f_high) or (f_low * f_high > 0):
                high *= 2.0
                f_high = beta_diff(high)
            expand_iter += 1

        # If still not bracketed, try symmetric contraction around estimate 0.2
        if not (np.isfinite(f_low) and np.isfinite(f_high) and f_low * f_high <= 0):
            # try scanning a grid and pick p with minimal absolute diff
            p_grid = np.logspace(-8, 3, 200)
            diffs = []
            for p in p_grid:
                d = beta_diff(p)
                diffs.append(np.abs(d) if np.isfinite(d) else np.inf)
            idx = int(np.argmin(diffs))
            p_opt = float(p_grid[idx])
        else:
            # Bisection
            tol = 1e-8
            max_iter = 80
            for _ in range(max_iter):
                mid = 0.5 * (low + high)
                f_mid = beta_diff(mid)
                if not np.isfinite(f_mid):
                    # shrink interval conservatively
                    high = 0.5 * (mid + high)
                    continue
                if abs(f_mid) < tol:
                    p_opt = mid
                    break
                if f_low * f_mid <= 0:
                    high = mid
                    f_high = f_mid
                else:
                    low = mid
                    f_low = f_mid
            else:
                # fallback to midpoint if loop exits without break
                p_opt = mid

    # Build the plasma object with the solved pressure parameter
    plasma = make_plasma(p_opt)

    # Plot pressure and F functions at the midplane.
    r_plot, z_plot = plasma.plotting_rz_arrays()

    psi_n_midplane = plasma.psi_norm(r_plot, 0)

    axis.set_xlabel(r"Radius [m]")
    axis.set_ylabel(r"Height [m]")
    axis.set_aspect("equal")

    psi_n = plasma.psi_norm(*np.meshgrid(r_plot, z_plot, indexing="ij"))

    # Plot contours of normalised poloidal flux.
    c = axis.contourf(
        r_plot, z_plot, psi_n.T, levels=np.linspace(0, 1.2, 13), cmap="plasma"
    )
    axis.contour(r_plot, z_plot, psi_n.T, colors="black", levels=np.linspace(0.1, 1, 10))
    fig.colorbar(c, ax=axis, label="Normalised Poloidal Flux")

    ax_equil_full = fig.add_subplot(2, 4, 1, aspect="equal")
    psi = plasma.psi(*np.meshgrid(r_plot, z_plot, indexing="ij"))
    # Mask out non-negative flux so only psi < 0 is plotted
    psi_masked = np.ma.masked_where(psi >= 0.0, psi)

    if psi_masked.mask.all():
        # Nothing negative to plot — show a message instead of an empty contour
        ax_equil_full.text(
            0.5,
            0.5,
            "No negative flux values to display",
            ha="center",
            va="center",
            transform=ax_equil_full.transAxes,
        )
    else:
        # build levels from the minimum negative value up to 0
        levels = np.linspace(float(psi_masked.min()), 0.0, 50)
        c_full = ax_equil_full.contourf(
            r_plot, z_plot, psi_masked.T, levels=levels, cmap="viridis"
        )
        fig.colorbar(c_full, ax=ax_equil_full, label="Poloidal Flux [Wb]")

    # keep existing normalized-psi contours (if desired)
    ax_equil_full.contour(
        r_plot, z_plot, psi_n.T, colors="black", levels=np.linspace(0.1, 1, 10)
    )

    ax_equil_full.set_xlabel(r"Radius [m]")
    ax_equil_full.set_ylabel(r"Height [m]")
    ax_equil_full.minorticks_on()
    # mark major radius and magnetic axis
    ax_equil_full.axvline(
        plasma.rmajor,
        color="black",
        linestyle="--",
        linewidth=1.2,
        label="Major radius $R_0$",
        zorder=150,
    )
    ax_equil_full.scatter(
        plasma.magnetic_axis[0],
        plasma.magnetic_axis[1],
        color="red",
        s=50,
        edgecolor="black",
        linewidth=0.6,
        zorder=200,
        label="Magnetic axis",
    )
    ax_equil_full.legend(loc="upper right", fontsize=8)
    ax_equil_full.grid(True, linestyle="--", alpha=0.5)

    ax_q_profile = fig.add_subplot(4, 2, 7)
    ax_q_profile.set_position([0.1, 0.075, 0.35, 0.15])  # [left, bottom, width, height]
    ax_q_profile.plot(np.linspace(0, 1, len(plasma.q_profile)), plasma.q_profile)
    ax_q_profile.set_xlabel(r"Normalised Radius $\rho$")
    ax_q_profile.set_ylabel("Safety Factor $q$")
    ax_q_profile.set_title("Equilibrium Safety Factor Profile")
    ax_q_profile.grid(True, linestyle="--", alpha=0.5)
    ax_q_profile.minorticks_on()
    ax_q_profile.set_xlim(0, 1)

    ax_equil_pressure = fig.add_subplot(4, 2, 5)
    ax_equil_pressure.set_position([
        0.1,
        0.3,
        0.35,
        0.15,
    ])  # [left, bottom, width, height]
    ax_equil_pressure.plot(
        r_plot,
        plasma.pressure_kpa(psi_n_midplane),
        color="black",
        label="Equilibrium Solved",
    )
    ax_equil_pressure.set_xlabel("Radius [m]")
    ax_equil_pressure.set_ylabel("Pressure [kPa]")
    ax_equil_pressure.set_title("Equilibrium Pressure Profile")
    ax_equil_pressure.grid(True, linestyle="--", alpha=0.5)
    ax_equil_pressure.minorticks_on()
    # overplot PROCESS pressure profile for comparison
    rho_profile = np.linspace(0, 1, n_plasma_profile_elements)
    r_profile = rminor * rho_profile + plasma.rmajor
    ax_equil_pressure.plot(
        r_profile,
        np.array(pres_plasma_thermal_total_profile) / 1e3,
        color="red",
        linestyle="--",
        label="PROCESS Profile",
    )
    # Mirror the plot in the x axis
    ax_equil_pressure.plot(
        2 * plasma.rmajor - r_profile,
        np.array(pres_plasma_thermal_total_profile) / 1e3,
        color="red",
        linestyle="--",
    )

    ax_equil_pressure.legend()

    # Show key 0D plasma parameters in an on-figure info box
    # compute relative and percentage differences (showing explicit + / -)
    if (
        np.isfinite(beta_thermal_poloidal_vol_avg)
        and abs(beta_thermal_poloidal_vol_avg) > 0
    ):
        rel_diff_bp = (
            plasma.beta_poloidal - beta_thermal_poloidal_vol_avg
        ) / beta_thermal_poloidal_vol_avg
    else:
        rel_diff_bp = np.nan
    pct_diff_bp = rel_diff_bp * 100 if np.isfinite(rel_diff_bp) else np.nan

    if (
        np.isfinite(beta_thermal_toroidal_vol_avg)
        and abs(beta_thermal_toroidal_vol_avg) > 0
    ):
        rel_diff_bt = (
            plasma.beta_toroidal - beta_thermal_toroidal_vol_avg
        ) / beta_thermal_toroidal_vol_avg
    else:
        rel_diff_bt = np.nan
    pct_diff_bt = rel_diff_bt * 100 if np.isfinite(rel_diff_bt) else np.nan

    # Use PROCESS reported total beta as the reference
    if np.isfinite(beta_thermal_vol_avg) and abs(beta_thermal_vol_avg) > 0:
        rel_diff_btot = (plasma.beta_total - beta_thermal_vol_avg) / beta_thermal_vol_avg
    else:
        rel_diff_btot = np.nan
    pct_diff_btot = rel_diff_btot * 100 if np.isfinite(rel_diff_btot) else np.nan

    if np.isfinite(beta_norm_thermal * 100) and abs(beta_norm_thermal) > 0:
        rel_diff_bnorm = (
            plasma.beta_normalised * 100 - beta_norm_thermal
        ) / beta_norm_thermal
    else:
        rel_diff_bnorm = np.nan
    pct_diff_bnorm = rel_diff_bnorm * 100 if np.isfinite(rel_diff_bnorm) else np.nan

    textstr = (
        f"Resolved pressure_parameter, $A$ = {p_opt:.3f}\n"
        f"Resolved kink safety factor = {plasma.kink_safety_factor:.3f}\n"
        f"Normalised circumference = {plasma.normalised_circumference:.3f}\n"
        f"Normalised volume = {plasma.normalised_volume:.3f}\n"
        # show beta poloidal and the signed relative & percentage difference vs PROCESS value
        f"Equilibria Beta poloidal = {plasma.beta_poloidal:.6f}\n"
        f"  Relative to PROCESS poloidal = {rel_diff_bp:+.3f}  ({pct_diff_bp:+.2f}%)\n"
        f"Equilibria Beta toroidal = {plasma.beta_toroidal:.6f}\n"
        f"  Relative to PROCESS toroidal = {rel_diff_bt:+.3f}  ({pct_diff_bt:+.2f}%)\n"
        f"Equilibria Beta total = {plasma.beta_total:.6f}\n"
        f"  Relative to PROCESS total = {rel_diff_btot:+.3f}  ({pct_diff_btot:+.2f}%)\n"
        f"Beta normalised = {plasma.beta_normalised:.6f}\n"
        f"  Relative to PROCESS normalised = {rel_diff_bnorm:+.3f}  ({pct_diff_bnorm:+.2f}%)\n"
        f"$\\Psi_0$ = {plasma.psi_0:.3f} Wb\n"
        f"$\\Psi_{{axis}}$ = {plasma.psi_axis:.3f} Wb\n"
        f"Magnetic axis = ({plasma.magnetic_axis[0]:.3f}, {plasma.magnetic_axis[1]:.3f}) m\n"
        f"Shafranov shift [m] = {plasma.dr_shafranov:.3f}"
    )

    # Place the info box in the top-left of the axis
    fig.text(
        0.35,
        0.85,
        textstr,
        fontsize=9,
        va="top",
        ha="left",
        bbox={"boxstyle": "round", "facecolor": "lightyellow", "alpha": 0.95},
    )

    # mark major radius and magnetic axis
    axis.axvline(
        plasma.rmajor,
        color="black",
        linestyle="--",
        linewidth=1.2,
        label="Major radius $R_0$",
        zorder=150,
    )
    axis.scatter(
        plasma.magnetic_axis[0],
        plasma.magnetic_axis[1],
        color="red",
        s=50,
        edgecolor="black",
        linewidth=0.6,
        zorder=200,
        label="Magnetic axis",
    )
    axis.legend(loc="upper right", fontsize=8)

    ax_poloidal_field = fig.add_subplot(4, 2, 6)
    ax_poloidal_field.set_position([
        0.6,
        0.3,
        0.35,
        0.15,
    ])  # [left, bottom, width, height]

    ax_poloidal_field.plot(
        r_plot,
        plasma.toroidal_current_density_ka_per_m2(r=r_plot),
        color="black",
        label="Equilibrium Solved",
    )
    ax_poloidal_field.set_xlabel("Radius [m]")
    ax_poloidal_field.grid(True, linestyle="--", alpha=0.5)
    ax_poloidal_field.minorticks_on()
    ax_poloidal_field.set_ylabel("Toroidal Current Density [kA/m$^2$]")
    ax_poloidal_field.set_title("Equilibrium Toroidal Current Density Profile")

    # Add title to the equilibrium analysis page
    fig.suptitle("Solov'ev Profiles Equilibrium Analysis", fontsize=16, y=0.95)
