def calculate_iter1996_nominal(dene20: float, bt: float, rmajor: float) -> float:
    """
    Calculate the nominal ITER-1996 L-H transition power threshold.

    :param dene20: Volume averaged electron density in units of 10^20 m^-3.
    :type dene20: float
    :param bt: Toroidal magnetic field [T]
    :type bt: float
    :param rmajor: Plasma major radius [m]
    :type rmajor: float
    :return: The ITER-1996 L-H transition power threshold [MW]
    :rtype: float

    :notes:

    :references:
        - T. Takizuka and International Atomic Energy Agency, Vienna (Austria),
        "Threshold power and energy confinement for ITER". 1996.

        - J. C. Wesley, “International Thermonuclear Experimental Reactor: Physics issues, capabilities and physics program plans,”
        Physics of Plasmas, vol. 4, no. 7, pp. 2642-2652, Jul. 1997,
        doi: https://doi.org/10.1063/1.872406.
    """
    return 0.45 * dene20**0.75 * bt * rmajor**2


def calculate_iter1996_upper(dene20: float, bt: float, rmajor: float) -> float:
    """
    Calculate the upper variant ITER-1996 L-H transition power threshold.

    :param dene20: Volume averaged electron density in units of 10^20 m^-3.
    :type dene20: float
    :param bt: Toroidal magnetic field [T]
    :type bt: float
    :param rmajor: Plasma major radius [m]
    :type rmajor: float
    :return: The ITER-1996 L-H transition power threshold [MW]
    :rtype: float

    :notes:

    :references:
        - T. Takizuka and International Atomic Energy Agency, Vienna (Austria),
        "Threshold power and energy confinement for ITER". 1996.

        - J. C. Wesley, “International Thermonuclear Experimental Reactor: Physics issues, capabilities and physics program plans,”
        Physics of Plasmas, vol. 4, no. 7, pp. 2642-2652, Jul. 1997,
        doi: https://doi.org/10.1063/1.872406.
    """
    return 0.3960502816 * dene20 * bt * rmajor**2.5


def calculate_iter1996_lower(dene20: float, bt: float, rmajor: float) -> float:
    """
    Calculate the lower variant ITER-1996 L-H transition power threshold.

    :param dene20: Volume averaged electron density in units of 10^20 m^-3.
    :type dene20: float
    :param bt: Toroidal magnetic field [T]
    :type bt: float
    :param rmajor: Plasma major radius [m]
    :type rmajor: float
    :return: The ITER-1996 L-H transition power threshold [MW]
    :rtype: float

    :notes:

    :references:
        - T. Takizuka and International Atomic Energy Agency, Vienna (Austria),
        "Threshold power and energy confinement for ITER". 1996.

        - J. C. Wesley, “International Thermonuclear Experimental Reactor: Physics issues, capabilities and physics program plans,”
        Physics of Plasmas, vol. 4, no. 7, pp. 2642-2652, Jul. 1997,
        doi: https://doi.org/10.1063/1.872406.
    """
    return 0.5112987149 * dene20**0.5 * bt * rmajor**1.5


def calculate_martin08_nominal(
    dnla20: float, bt: float, a_plasma_surface: float, m_ions_total_amu: float
) -> float:
    """
        Calculate the nominal Martin L-H transition power threshold.

        :param dnla20: Line averaged electron density in units of 10^20 m^-3.
        :type dnla20: float
        :param bt: Toroidal magnetic field [T]
        :type bt: float
        :param a_plasma_surface: Plasma surface area [m^2]
        :type a_plasma_surface: float
        :param m_ions_total_amu: Total ion mass in atomic mass units [amu]
        :type m_ions_total_amu: float
        :return: The Martin L-H transition power threshold [MW]
        :rtype: float

        :notes:

        :references:
            - Y. R. Martin, T. Takizuka, and the I. C. H-mode. T. D. Group, “Power requirement for accessing the H-mode in ITER,”
            Journal of Physics: Conference Series, vol. 123, p. 012033, Jul. 2008,
            doi: https://doi.org/10.1088/1742-6596/123/1/012033.
    ‌
    """
    return (
        0.0488
        * dnla20**0.717
        * bt**0.803
        * a_plasma_surface**0.941
        * (2.0 / m_ions_total_amu)
    )


def calculate_martin08_upper(
    dnla20: float, bt: float, a_plasma_surface: float, m_ions_total_amu: float
) -> float:
    """
        Calculate the upper Martin L-H transition power threshold.

        :param dnla20: Line averaged electron density in units of 10^20 m^-3.
        :type dnla20: float
        :param bt: Toroidal magnetic field [T]
        :type bt: float
        :param a_plasma_surface: Plasma surface area [m^2]
        :type a_plasma_surface: float
        :param m_ions_total_amu: Total ion mass in atomic mass units [amu]
        :type m_ions_total_amu: float
        :return: The Martin L-H transition power threshold [MW]
        :rtype: float

        :notes:

        :references:
            - Y. R. Martin, T. Takizuka, and the I. C. H-mode. T. D. Group, “Power requirement for accessing the H-mode in ITER,”
            Journal of Physics: Conference Series, vol. 123, p. 012033, Jul. 2008,
            doi: https://doi.org/10.1088/1742-6596/123/1/012033.
    ‌
    """
    return (
        0.05166240355
        * dnla20**0.752
        * bt**0.835
        * a_plasma_surface**0.96
        * (2.0 / m_ions_total_amu)
    )


def calculate_martin08_lower(
    dnla20: float, bt: float, a_plasma_surface: float, m_ions_total_amu: float
) -> float:
    """
        Calculate the lower Martin L-H transition power threshold.

        :param dnla20: Line averaged electron density in units of 10^20 m^-3.
        :type dnla20: float
        :param bt: Toroidal magnetic field [T]
        :type bt: float
        :param a_plasma_surface: Plasma surface area [m^2]
        :type a_plasma_surface: float
        :param m_ions_total_amu: Total ion mass in atomic mass units [amu]
        :type m_ions_total_amu: float
        :return: The Martin L-H transition power threshold [MW]
        :rtype: float

        :notes:

        :references:
            - Y. R. Martin, T. Takizuka, and the I. C. H-mode. T. D. Group, “Power requirement for accessing the H-mode in ITER,”
            Journal of Physics: Conference Series, vol. 123, p. 012033, Jul. 2008,
            doi: https://doi.org/10.1088/1742-6596/123/1/012033.
    ‌
    """
    return (
        0.04609619059
        * dnla20**0.682
        * bt**0.771
        * a_plasma_surface**0.922
        * (2.0 / m_ions_total_amu)
    )
