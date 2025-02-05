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
