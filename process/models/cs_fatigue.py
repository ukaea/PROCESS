"""Module containing CS fatigue routines"""

import numpy as np
from numba import njit

from process.core import constants
from process.core import process_output as op
from process.core.model import Model


class CSFatigue(Model):
    """Calculate CS fatigue model parameters"""

    def __init__(self):
        self.outfile = constants.NOUT

    def output(self):
        """Output CS fatigue model parameters to the output file"""
        op.osubhd(self.outfile, "CS Cyclical Stress :")
        if self.data.physics.f_c_plasma_inductive > 0.0e-4:
            op.ovarre(
                self.outfile,
                "Residual hoop stress in CS Steel (Pa)",
                "(stress_hoop_cs_residual)",
                self.data.cs_fatigue.stress_hoop_cs_residual,
            )
            op.ovarre(
                self.outfile,
                "Minimum burn time (s)",
                "(t_burn_min)",
                self.data.constraints.t_burn_min,
            )
            op.ovarre(
                self.outfile,
                "Initial vertical crack size (m)",
                "(dz_cs_turn_crack_initial)",
                self.data.cs_fatigue.dz_cs_turn_crack_initial,
            )
            op.ovarre(
                self.outfile,
                "Initial radial crack size (m)",
                "(dr_cs_turn_crack_initial)",
                self.data.cs_fatigue.dr_cs_turn_crack_initial,
            )

            op.ovarre(
                self.outfile,
                "Allowable number of cycles till CS fracture",
                "(n_cycle)",
                self.data.cs_fatigue.n_cycle,
                "OP ",
            )
            op.ovarre(
                self.outfile,
                "Minimum number of cycles required till CS fracture",
                "(n_cycle_min)",
                self.data.cs_fatigue.n_cycle_min,
                "OP ",
            )

    def run(self):
        """CSFatigue model doesn't need to be run"""

    def ncycle(
        self,
        max_hoop_stress: float,
        residual_stress: float,
        dz_cs_turn_crack_initial: float,
        dz_cs_turn_conduit: float,
        dr_cs_turn_conduit: float,
    ) -> tuple[float, float]:
        """

        Parameters
        ----------
        max_hoop_stress :

        residual_stress :

        dz_cs_turn_crack_initial :

        dz_cs_turn_conduit :

        dr_cs_turn_conduit :

        """
        # Default Parameters for SS 316LN from
        # X. Sarasola et al, IEEE Transactions on Applied Superconductivity,
        # vol. 30, no. 4, pp. 1-5, June 2020, Art no. 4200705

        n = -self.data.cs_fatigue.paris_power_law * (
            self.data.cs_fatigue.walker_coefficient - 1.0e0
        )

        # Set units to MPa
        max_hoop_stress_MPa = max_hoop_stress / 1.0e6
        residual_stress_MPa = residual_stress / 1.0e6

        # Set initial crack size
        dr_cs_turn_crack_initial = 3.0e0 * dz_cs_turn_crack_initial
        a = dz_cs_turn_crack_initial
        c = dr_cs_turn_crack_initial

        # Cyclic element of stress
        hoop_stress_MPa = max_hoop_stress_MPa

        # Mean stress ratio
        # Fatigue Stress Assessment in Fusion Magnet Components
        # J. Lorenzo, X. Sarasola, M. Mantsinen
        r = residual_stress_MPa / (max_hoop_stress_MPa + residual_stress_MPa)

        # Calculated constant for a given stress ratio using Walker equation
        # https://en.wikipedia.org/wiki/Crack_growth_equation#Walker_equation
        cr = self.data.cs_fatigue.paris_coefficient_cs_turn / (1.0e0 - r) ** n

        # select given increase in crack area
        delta = 1.0e-4

        # Initialise number of cycles
        n_pulse = 0.0
        k_max = 0.0

        # factor 2 taken as safety factors in the crack sizes
        # Default CS steel undergoes fast fracture when SIF > 200 MPa,
        # under a safety factor 1.5 we use 133MPa
        pi_2_arr = np.array([np.pi / 2.0e0, 0])
        while (
            (a <= dz_cs_turn_conduit / self.data.cs_fatigue.sf_vertical_crack)
            and (c <= dr_cs_turn_conduit / self.data.cs_fatigue.sf_radial_crack)
            and (
                k_max
                <= self.data.cs_fatigue.fracture_toughness
                / self.data.cs_fatigue.sf_fast_fracture
            )
        ):
            # find SIF max from SIF_a and SIF_c
            k_a, k_c = self.surface_stress_intensity_factor(
                hoop_stress_MPa,
                dz_cs_turn_conduit,
                dr_cs_turn_conduit,
                a,
                c,
                pi_2_arr,
            )
            k_max = max(k_a, k_c)

            # run euler_method and find number of cycles needed to give crack increase
            delta_n = delta / (cr * (k_max**self.data.cs_fatigue.paris_power_law))

            # update a and c, N (+= doesn't work for fortran (?) reasons)
            a += delta * (k_a / k_max) ** self.data.cs_fatigue.paris_power_law
            c += delta * (k_c / k_max) ** self.data.cs_fatigue.paris_power_law
            n_pulse += delta_n

        # two pulses - ramp to Vsmax and ramp down per cycle
        return n_pulse / 2.0e0, dr_cs_turn_crack_initial

    @staticmethod
    @njit(cache=True)
    def embedded_stress_intensity_factor(hoop_stress, t, w, a, c, phi):
        """
        Assumes an embedded elliptical effect geometry

        geometric quantities
        hoop_stress - change in hoop stress over cycle
        t - plate thickness
        w - plate width
        a - crack depth (t -direction)
        c - crack length (w - direction)
        Ref: J. C. Newman, I. S. Raju "Stress-Intensity Factor Equations for Cracks in
        Three-Dimensional Finite Bodies Subjected to Tension and Bending Loads"
        https://core.ac.uk/download/pdf/42849129.pdf
        Ref: C. Jong, Magnet Structural Design
        Criteria Part 1: Main Structural Components and Welds 2012
        """
        # reuse of calc
        a_c = a / c
        a_t = a / t
        cos_phi = np.cos(phi)
        cos_phi_2 = cos_phi**2.0e0
        sin_phi_2 = np.sin(phi) ** 2.0e0

        if a <= c:
            q = 1.0e0 + 1.464e0 * a_c**1.65e0
            m1 = 1.0e0
            f_phi = (a_c**2.0e0 * cos_phi_2 + sin_phi_2) ** 0.25e0
        else:  # elif a > c:
            c_a = c / a
            q = 1.0e0 + 1.464e0 * c_a**1.65e0
            m1 = np.sqrt(c_a)
            f_phi = (c_a**2.0e0 * sin_phi_2 + cos_phi_2) ** 0.25e0

        # compute the unitless geometric correction
        # compute the stress intensity factor
        return (
            hoop_stress
            * (  # f
                (
                    m1
                    + (0.05e0 / (0.11e0 + a_c**1.5e0)) * a_t**2.0e0  # m2 *a_t^2
                    + (0.29e0 / (0.23e0 + a_c**1.5e0)) * a_t**4.0e0  # m3 *a_t^4
                )
                * (  # g
                    1.0e0
                    - (
                        a_t**4.0e0
                        * np.sqrt(2.6e0 - (2.0e0 * a_t))
                        / (1.0e0 + 4.0e0 * a_c)
                    )
                    * abs(cos_phi)
                )
                * f_phi
                * np.sqrt(  # f_w
                    1.0e0 / np.cos(np.sqrt(a_t) * np.pi * c / (2.0e0 * w))
                )
            )
            * np.sqrt(np.pi * a / q)
        )

    @staticmethod
    @njit(cache=True)
    def surface_stress_intensity_factor(hoop_stress, t, w, a, c, phi):
        """
        Calculate the stress intensity factor (K) for a semi-elliptical surface crack in
        a finite plate under tension and bending loads.

        geometric quantities
        hoop_stress - change in hoop stress over cycle
        t - plate thickness
        w - plate width
        a - crack depth (t -direction)
        c - crack length (w - direction)


        Ref: C. Jong, Magnet Structural Design
        Criteria Part 1: Main Structural Components and Welds 2012

        Notes
        -----
        - The semi-elliptical surface crack example can be seen in Figure 2(b) of
        Reference [1]

        References
        ----------
        [1] Newman, Isadore & Raju, Ivatury. (1984). Stress Intensity Factor Equations
        for Cracks in Three-Dimensional Finite Bodies Subjected to Tension and Bending
        Loads. NASA Technical Memorandum. 85.
        https://ntrs.nasa.gov/citations/19840015857
        """
        bending_stress = 0.0e0  # * 3.0 * M / (w*d**2.0)

        # Depth of crack is less than or equal to the half-length of the crack
        if a <= c:
            # Equation 24 of Reference [1]
            g_21 = -1.22e0 - 0.12e0 * (a / c)

            # Equation 25 of Reference [1]
            g_22 = 0.55e0 - 1.05e0 * (a / c) ** 0.75e0 + 0.47e0 * (a / c) ** 1.5e0

            # Equation 2(a) of Reference [1]
            q = 1.0e0 + 1.464e0 * (a / c) ** 1.65e0

            # Equation 16 of Reference [1]
            m1 = 1.13e0 - 0.09e0 * (a / c)

            # Equation 17 of Reference [1]
            m2 = -0.54e0 + 0.89e0 / (0.2e0 + (a / c))

            # Equation 18 of Reference [1]
            m3 = 0.5e0 - 1.0e0 / (0.65e0 + (a / c)) + 14.0e0 * (1 - (a / c)) ** 24.0e0

            # Equation 19 of Reference [1]
            g = (
                1.0e0
                + (0.1e0 + 0.35e0 * (a / t) ** 2.0e0) * (1.0e0 - np.sin(phi)) ** 2.0e0
            )

            # Equation 10 of Reference [1]
            f_phi = (
                (a / c) ** 2.0e0 * (np.cos(phi) ** 2.0e0) + np.sin(phi) ** 2.0e0
            ) ** 0.25e0

            # Equation 21 of Reference [1]
            p = 0.2e0 + (a / c) + 0.6e0 * (a / t)

            # Equation 22 of Reference [1]
            H1 = 1.0e0 - 0.34e0 * (a / t) - 0.11e0 * a * a / (c * t)

            # Equation 23 of Reference [1]
            H2 = 1.0 + (g_21) * (a / t) + (g_22) * (a / t) ** 2.0

        # Depth of crack is greater than the half-length of the crack (a > c)
        else:
            # Equation 33 of Reference [1]
            g_11 = -0.04 - 0.41 * (c / a)

            # Equation 34 of Reference [1]
            g_12 = 0.55 - 1.93 * (c / a) ** 0.75 + 1.38 * (c / a) ** 1.5

            # Equation 35 of Reference [1]
            g_21 = -2.11 + (0.77 * (c / a))

            # Equation 36 of Reference [1]
            g_22 = 0.55 - 0.72 * (c / a) ** 0.75 + 0.14 * (c / a) ** 1.5

            # Equation 2(b) of Reference [1]
            q = 1.0e0 + 1.464e0 * (c / a) ** 1.65e0

            # Equation 26 of Reference [1]
            m1 = np.sqrt(c / a) * (1.0e0 + 0.04e0 * (c / a))

            # Equation 27 of Reference [1]
            m2 = 0.2e0 * (c / a) ** 4.0

            # Equation 28 of Reference [1]
            m3 = -0.11e0 * (c / a) ** 4.0

            # Equation 29 of Reference [1]
            g = (
                1.0e0
                + (0.1e0 + 0.35e0 * (c / a) * (a / t) ** 2.0)
                * (1.0e0 - np.sin(phi)) ** 2.0e0
            )

            # Equation 13 of Reference [1]
            f_phi = (
                (c / a) ** 2.0e0 * np.sin(phi) ** 2.0e0 + (np.cos(phi) ** 2.0e0)
            ) ** 0.25e0

            # Equation 30 of Reference [1]
            p = 0.2e0 + (c / a) + 0.6e0 * (a / t)

            # Equation 31 of Reference [1]
            H1 = 1.0e0 + (g_11) * (a / t) + (g_12) * (a / t) ** 2.0

            # Equation 32 of Reference [1]
            H2 = 1.0e0 + (g_21) * (a / t) + (g_22) * (a / t) ** 2.0

        # Equation 11 of Reference [1]
        # Finite width correction factor
        f_w = np.sqrt(1.0e0 / np.cos(np.sqrt(a / t) * np.pi * c / (2.0e0 * w)))

        # Equation 15 of Reference [1]
        # Boundary-correction factor for surface crack in a plate under tension
        f_s = (m1 + m2 * (a / t) ** 2.0 + m3 * (a / t) ** 4.0e0) * g * f_phi * f_w

        # Equation 1(c) of Reference [1]
        # Bending multiplier for surface crack in a plate
        h_s = H1 + (H2 - H1) * np.sin(phi) ** p

        # Equation 1(a) of Reference [1]
        # Stress intensity factor for a semi-elliptical surface crack in a finite plate
        return (hoop_stress + (h_s) * bending_stress) * f_s * np.sqrt(np.pi * a / q)
