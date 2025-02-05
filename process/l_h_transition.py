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


def calculate_snipes2000_nominal(
    dnla20: float, bt: float, rmajor: float, rminor: float, m_ions_total_amu: float
) -> float:
    """
        Calculate the nominal Snipes 2000 L-H transition power threshold.

        :param dnla20: Line averaged electron density in units of 10^20 m^-3.
        :type dnla20: float
        :param bt: Toroidal magnetic field [T]
        :type bt: float
        :param rmajor: Plasma major radius [m]
        :type rmajor: float
        :param rminor: Plasma minor radius [m]
        :type rminor: float
        :param m_ions_total_amu: Total ion mass in atomic mass units [amu]
        :type m_ions_total_amu: float
        :return: The Snipes 2000 L-H transition power threshold [MW]
        :rtype: float

        :notes:

        :references:
            - J. A. Snipes and the I. H-mode. T. Group, “Latest results on the H-mode threshold using the international H-mode threshold database,”
            Plasma Physics and Controlled Fusion, vol. 42, no. 5A, pp. A299-A308, May 2000,
            doi: https://doi.org/10.1088/0741-3335/42/5a/336.
    ‌
    """
    return (
        1.42
        * dnla20**0.58
        * bt**0.82
        * rmajor
        * rminor**0.81
        * (2.0 / m_ions_total_amu)
    )


def calculate_snipes2000_upper(
    dnla20: float, bt: float, rmajor: float, rminor: float, m_ions_total_amu: float
) -> float:
    """
        Calculate the upper Snipes 2000 L-H transition power threshold.

        :param dnla20: Line averaged electron density in units of 10^20 m^-3.
        :type dnla20: float
        :param bt: Toroidal magnetic field [T]
        :type bt: float
        :param rmajor: Plasma major radius [m]
        :type rmajor: float
        :param rminor: Plasma minor radius [m]
        :type rminor: float
        :param m_ions_total_amu: Total ion mass in atomic mass units [amu]
        :type m_ions_total_amu: float
        :return: The Snipes 2000 L-H transition power threshold [MW]
        :rtype: float

        :notes:

        :references:
            - J. A. Snipes and the I. H-mode. T. Group, “Latest results on the H-mode threshold using the international H-mode threshold database,”
            Plasma Physics and Controlled Fusion, vol. 42, no. 5A, pp. A299-A308, May 2000,
            doi: https://doi.org/10.1088/0741-3335/42/5a/336.
    ‌
    """
    return (
        1.547
        * dnla20**0.615
        * bt**0.851
        * rmajor**1.089
        * rminor**0.876
        * (2.0 / m_ions_total_amu)
    )


def calculate_snipes2000_lower(
    dnla20: float, bt: float, rmajor: float, rminor: float, m_ions_total_amu: float
) -> float:
    """
        Calculate the lower Snipes 2000 L-H transition power threshold.

        :param dnla20: Line averaged electron density in units of 10^20 m^-3.
        :type dnla20: float
        :param bt: Toroidal magnetic field [T]
        :type bt: float
        :param rmajor: Plasma major radius [m]
        :type rmajor: float
        :param rminor: Plasma minor radius [m]
        :type rminor: float
        :param m_ions_total_amu: Total ion mass in atomic mass units [amu]
        :type m_ions_total_amu: float
        :return: The Snipes 2000 L-H transition power threshold [MW]
        :rtype: float

        :notes:

        :references:
            - J. A. Snipes and the I. H-mode. T. Group, “Latest results on the H-mode threshold using the international H-mode threshold database,”
            Plasma Physics and Controlled Fusion, vol. 42, no. 5A, pp. A299-A308, May 2000,
            doi: https://doi.org/10.1088/0741-3335/42/5a/336.
    ‌
    """
    return (
        1.293
        * dnla20**0.545
        * bt**0.789
        * rmajor**0.911
        * rminor**0.744
        * (2.0 / m_ions_total_amu)
    )


def calculate_snipes2000_closed_divertor_nominal(
    dnla20: float, bt: float, rmajor: float, m_ions_total_amu: float
) -> float:
    """
        Calculate the nominal Snipes 2000 Closed Divertor L-H transition power threshold with CD factor.

        :param dnla20: Line averaged electron density in units of 10^20 m^-3.
        :type dnla20: float
        :param bt: Toroidal magnetic field [T]
        :type bt: float
        :param rmajor: Plasma major radius [m]
        :type rmajor: float
        :param m_ions_total_amu: Total ion mass in atomic mass units [amu]
        :type m_ions_total_amu: float
        :return: The Snipes 2000 L-H transition power threshold with CD factor [MW]
        :rtype: float

        :notes:

        :references:
            - J. A. Snipes and the I. H-mode. T. Group, “Latest results on the H-mode threshold using the international H-mode threshold database,”
            Plasma Physics and Controlled Fusion, vol. 42, no. 5A, pp. A299-A308, May 2000,
            doi: https://doi.org/10.1088/0741-3335/42/5a/336.
    ‌
    """
    return 0.8 * dnla20**0.5 * bt**0.53 * rmajor**1.51 * (2.0 / m_ions_total_amu)


def calculate_snipes2000_closed_divertor_upper(
    dnla20: float, bt: float, rmajor: float, m_ions_total_amu: float
) -> float:
    """
        Calculate the upper Snipes 2000 Closed Divertor L-H transition power threshold with CD factor.

        :param dnla20: Line averaged electron density in units of 10^20 m^-3.
        :type dnla20: float
        :param bt: Toroidal magnetic field [T]
        :type bt: float
        :param rmajor: Plasma major radius [m]
        :type rmajor: float
        :param m_ions_total_amu: Total ion mass in atomic mass units [amu]
        :type m_ions_total_amu: float
        :return: The Snipes 2000 L-H transition power threshold with CD factor [MW]
        :rtype: float

        :notes:

        :references:
            - J. A. Snipes and the I. H-mode. T. Group, “Latest results on the H-mode threshold using the international H-mode threshold database,”
            Plasma Physics and Controlled Fusion, vol. 42, no. 5A, pp. A299-A308, May 2000,
            doi: https://doi.org/10.1088/0741-3335/42/5a/336.
    ‌
    """
    return 0.867 * dnla20**0.561 * bt**0.588 * rmajor**1.587 * (2.0 / m_ions_total_amu)


def calculate_snipes2000_closed_divertor_lower(
    dnla20: float, bt: float, rmajor: float, m_ions_total_amu: float
) -> float:
    """
        Calculate the lower Snipes 2000 Closed Divertor L-H transition power threshold with CD factor.

        :param dnla20: Line averaged electron density in units of 10^20 m^-3.
        :type dnla20: float
        :param bt: Toroidal magnetic field [T]
        :type bt: float
        :param rmajor: Plasma major radius [m]
        :type rmajor: float
        :param m_ions_total_amu: Total ion mass in atomic mass units [amu]
        :type m_ions_total_amu: float
        :return: The Snipes 2000 L-H transition power threshold with CD factor [MW]
        :rtype: float

        :notes:

        :references:
            - J. A. Snipes and the I. H-mode. T. Group, “Latest results on the H-mode threshold using the international H-mode threshold database,”
            Plasma Physics and Controlled Fusion, vol. 42, no. 5A, pp. A299-A308, May 2000,
            doi: https://doi.org/10.1088/0741-3335/42/5a/336.
    ‌
    """
    return 0.733 * dnla20**0.439 * bt**0.472 * rmajor**1.433 * (2.0 / m_ions_total_amu)


def calculate_hubbard2012_nominal(plasma_current: float, dnla20: float) -> float:
    """
        Calculate the nominal Hubbard 2012 L-I transition power threshold.

        :param plasma_current: Plasma current [A]
        :type plasma_current: float
        :param dnla20: Line averaged electron density in units of 10^20 m^-3.
        :type dnla20: float
        :return: The Hubbard 2012 L-I transition power threshold [MW]
        :rtype: float

        :notes:

        :references:
            - A. E. Hubbard et al., “Threshold conditions for transitions to I-mode and H-mode with unfavourable ion grad B drift direction,”
            Nuclear Fusion, vol. 52, no. 11, pp. 114009-114009, Oct. 2012,
            doi: https://doi.org/10.1088/0029-5515/52/11/114009.
    ‌
    """
    return 2.11 * (plasma_current / 1e6) ** 0.94 * dnla20**0.65


def calculate_hubbard2012_upper(plasma_current: float, dnla20: float) -> float:
    """
        Calculate the upper Hubbard 2012 L-I transition power threshold.

        :param plasma_current: Plasma current [A]
        :type plasma_current: float
        :param dnla20: Line averaged electron density in units of 10^20 m^-3.
        :type dnla20: float
        :return: The Hubbard 2012 L-I transition power threshold [MW]
        :rtype: float

        :notes:

        :references:
            - A. E. Hubbard et al., “Threshold conditions for transitions to I-mode and H-mode with unfavourable ion grad B drift direction,”
            Nuclear Fusion, vol. 52, no. 11, pp. 114009-114009, Oct. 2012,
            doi: https://doi.org/10.1088/0029-5515/52/11/114009.
    ‌
    """
    return 2.11 * (plasma_current / 1e6) ** 1.18 * dnla20**0.83


def calculate_hubbard2012_lower(plasma_current: float, dnla20: float) -> float:
    """
        Calculate the lower Hubbard 2012 L-I transition power threshold.

        :param plasma_current: Plasma current [A]
        :type plasma_current: float
        :param dnla20: Line averaged electron density in units of 10^20 m^-3.
        :type dnla20: float
        :return: The Hubbard 2012 L-I transition power threshold [MW]
        :rtype: float

        :notes:

        :references:
            - A. E. Hubbard et al., “Threshold conditions for transitions to I-mode and H-mode with unfavourable ion grad B drift direction,”
            Nuclear Fusion, vol. 52, no. 11, pp. 114009-114009, Oct. 2012,
            doi: https://doi.org/10.1088/0029-5515/52/11/114009.
    ‌
    """
    return 2.11 * (plasma_current / 1e6) ** 0.7 * dnla20**0.47
