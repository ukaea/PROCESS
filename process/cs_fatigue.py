import numpy as np
from numba import njit

import process.data_structure as data_structure
from process import constants


class CsFatigue:
    def __init__(self):
        self.outfile = constants.NOUT

    def ncycle(
        self,
        max_hoop_stress,
        residual_stress,
        t_crack_vertical,
        dz_cs_turn_conduit,
        dr_cs_turn_conduit,
    ):
        """

        Parameters
        ----------
        max_hoop_stress :

        residual_stress :

        t_crack_vertical :

        dz_cs_turn_conduit :

        dr_cs_turn_conduit :

        """
        # Default Parameters for SS 316LN from
        # X. Sarasola et al, IEEE Transactions on Applied Superconductivity,
        # vol. 30, no. 4, pp. 1-5, June 2020, Art no. 4200705

        n = -data_structure.cs_fatigue_variables.paris_power_law * (
            data_structure.cs_fatigue_variables.walker_coefficient - 1.0e0
        )

        # Set units to MPa
        max_hoop_stress_MPa = max_hoop_stress / 1.0e6
        residual_stress_MPa = residual_stress / 1.0e6

        # Set intial crack size
        t_crack_radial = 3.0e0 * t_crack_vertical
        a = t_crack_vertical
        c = t_crack_radial

        # Cyclic element of stress
        hoop_stress_MPa = max_hoop_stress_MPa

        # Mean stress ratio
        # Fatigue Stress Assessment in Fusion Magnet Components
        # J. Lorenzo, X. Sarasola, M. Mantsinen
        r = residual_stress_MPa / (max_hoop_stress_MPa + residual_stress_MPa)

        # Calculated constant for a given stress ratio using Walker equation
        # https://en.wikipedia.org/wiki/Crack_growth_equation#Walker_equation
        cr = data_structure.cs_fatigue_variables.paris_coefficient / (1.0e0 - r) ** n

        # select given increase in crack area
        delta = 1.0e-4

        # Initialise number of cycles
        n_pulse = 0.0
        k_max = 0.0

        # factor 2 taken as saftey factors in the crack sizes
        # Default CS steel undergoes fast fracture when SIF > 200 MPa, under a saftey factor 1.5 we use 133MPa
        pi_2_arr = np.array([np.pi / 2.0e0, 0])
        while (
            (
                a
                <= dz_cs_turn_conduit
                / data_structure.cs_fatigue_variables.sf_vertical_crack
            )
            and (
                c
                <= dr_cs_turn_conduit
                / data_structure.cs_fatigue_variables.sf_radial_crack
            )
            and (
                k_max
                <= data_structure.cs_fatigue_variables.fracture_toughness
                / data_structure.cs_fatigue_variables.sf_fast_fracture
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
            delta_n = delta / (
                cr * (k_max**data_structure.cs_fatigue_variables.paris_power_law)
            )

            # update a and c, N (+= doesnt work for fortran (?) reasons)
            a = (
                a
                + delta
                * (k_a / k_max) ** data_structure.cs_fatigue_variables.paris_power_law
            )
            c = (
                c
                + delta
                * (k_c / k_max) ** data_structure.cs_fatigue_variables.paris_power_law
            )
            n_pulse = n_pulse + delta_n

        # two pulses - ramp to Vsmax and ramp down per cycle
        return n_pulse / 2.0e0, t_crack_radial

    @staticmethod
    @njit(cache=True)
    def embedded_stress_intensity_factor(hoop_stress, t, w, a, c, phi):
        # ! Assumes an embedded elliptical efect geometry
        # ! geometric quantities
        # ! hoop_stress - change in hoop stress over cycle
        # ! t - plate thickness
        # ! w - plate width
        # ! a - crack depth (t -direction)
        # ! c - crack length (w - direction)
        # ! Ref: J. C. Newman, I. S. Raju "Stress-Intensity Factor Equations for Cracks in
        # ! Three-Dimensional Finite Bodies Subjected to Tension and Bending Loads"
        # ! https://core.ac.uk/download/pdf/42849129.pdf
        # ! Ref: C. Jong, Magnet Structural Design
        # ! Criteria Part 1: Main Structural Components and Welds 2012

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
        # ! Assumes an surface semi elliptical defect geometry
        # ! geometric quantities
        # ! hoop_stress - change in hoop stress over cycle
        # ! t - plate thickness
        # ! w - plate width
        # ! a - crack depth (t -direction)
        # ! c - crack length (w - direction)
        # ! Ref: J. C. Newman, I. S. Raju "Stress-Intensity Factor Equations for Cracks in
        # ! Three-Dimensional Finite Bodies Subjected to Tension and Bending Loads"
        # ! https://core.ac.uk/download/pdf/42849129.pdf
        # ! Ref: C. Jong, Magnet Structural Design
        # ! Criteria Part 1: Main Structural Components and Welds 2012

        bending_stress = 0.0e0  # * 3.0 * M / (w*d**2.0)

        # reuse of calc
        a_t = a / t
        a_t_2 = a_t**2.0e0
        sin_phi = np.sin(phi)
        cos_phi_2 = np.cos(phi) ** 2.0e0

        if a <= c:
            # reuse of calc
            a_c = a / c

            q = 1.0e0 + 1.464e0 * a_c**1.65e0
            m1 = 1.13e0 - 0.09e0 * a_c
            m2 = -0.54e0 + 0.89e0 / (0.2e0 + a_c)
            m3 = 0.5e0 - 1.0e0 / (0.65e0 + a_c) + 14.0e0 * (1 - a_c) ** 24.0e0
            g = 1.0e0 + (0.1e0 + 0.35e0 * a_t**2.0e0) * (1.0e0 - sin_phi) ** 2.0e0
            f_phi = (a_c**2.0e0 * cos_phi_2 + sin_phi**2.0e0) ** 0.25e0
            p = 0.2e0 + a_c + 0.6e0 * a_t
            H1 = 1.0e0 - 0.34e0 * a_t - 0.11e0 * a * a / (c * t)
            H2 = (
                1.0
                + (-1.22e0 - 0.12e0 * a_c) * a_t  # G21 * a / t
                + (0.55e0 - 1.05e0 * a_c**0.75e0 + 0.47e0 * a_c**1.5e0)  # G22
                * a_t_2
            )
        else:  # elif a > c:
            c_a = c / a
            c_a_4 = c_a**4.0e0

            q = 1.0e0 + 1.464e0 * c_a**1.65e0
            m1 = np.sqrt(c_a) * (1.0e0 + 0.04e0 * c_a)

            m2 = 0.2e0 * c_a_4
            m3 = -0.11e0 * c_a_4

            g = 1.0e0 + (0.1e0 + 0.35e0 * c_a * a_t_2) * (1.0e0 - sin_phi) ** 2.0e0
            f_phi = (c_a**2.0e0 * sin_phi**2.0e0 + cos_phi_2) ** 0.25e0
            p = 0.2e0 + c_a + 0.6e0 * a_t
            H1 = (
                1.0e0
                + (-0.04e0 - 0.41e0 * c_a) * a_t  # G11 * a / t
                + (0.55e0 - 1.93e0 * c_a**0.75e0 + 1.38e0 * c_a**1.5e0)  # G12
                * a_t_2
            )
            H2 = (
                1.0e0
                + (-2.11e0 + 0.77e0 * c_a) * a_t  # G21 * a / t
                + (0.55e0 - 0.72e0 * c_a**0.75e0 + 0.14e0 * c_a * 1.5e0) * a_t_2  # G22
            )

        # compute the unitless geometric correction
        # compute the stress intensity factor
        return (
            (  # hoop_stress + Hs * bending_stress
                hoop_stress + (H1 + (H2 - H1) * sin_phi**p) * bending_stress
            )
            * (  # f
                (m1 + m2 * a_t_2 + m3 * a_t**4.0e0)
                * g
                * f_phi
                * np.sqrt(  # f_w
                    1.0e0 / np.cos(np.sqrt(a_t) * np.pi * c / (2.0e0 * w))
                )
            )
            * np.sqrt(np.pi * a / q)
        )
