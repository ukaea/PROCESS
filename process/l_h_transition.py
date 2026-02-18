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

        - J. C. Wesley, “International Thermonuclear Experimental Reactor: Physics issues, capabilities and physics program plans,”
        Physics of Plasmas, vol. 4, no. 7, pp. 2642-2652, Jul. 1997,
        doi: https://doi.org/10.1063/1.872406.
    """
    return 0.45 * dnla20**0.75 * b_plasma_toroidal_on_axis * rmajor**2


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

        - J. C. Wesley, “International Thermonuclear Experimental Reactor: Physics issues, capabilities and physics program plans,”
        Physics of Plasmas, vol. 4, no. 7, pp. 2642-2652, Jul. 1997,
        doi: https://doi.org/10.1063/1.872406.
    """
    return 0.3960502816 * dnla20 * b_plasma_toroidal_on_axis * rmajor**2.5


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

        - J. C. Wesley, “International Thermonuclear Experimental Reactor: Physics issues, capabilities and physics program plans,”
        Physics of Plasmas, vol. 4, no. 7, pp. 2642-2652, Jul. 1997,
        doi: https://doi.org/10.1063/1.872406.
    """
    return 0.5112987149 * dnla20**0.5 * b_plasma_toroidal_on_axis * rmajor**1.5


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
        - J. A. Snipes and the ITER H-mode Threshold Database Working Group, "An Analysis of the H-mode Threshold in ITER,"
        Controlled Fusion and Plasma Physics, 24th EPS Conference,
        Berchtesgaden, June 9th-13th 1997, vol.21A, part III, p.961.
        url:https://library.ipp.mpg.de/EPS_24_Vol3_1997.pdf.
        *This is a conference poster*
    """
    return 0.65 * dnla20**0.93 * b_plasma_toroidal_on_axis**0.86 * rmajor**2.15


def calculate_snipes1997_kappa(
    dnla20: float, b_plasma_toroidal_on_axis: float, rmajor: float, kappa: float
) -> float:
    """Calculate the Snipes 1997 ITER L-H transition power threshold with kappa factor.

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
        - J. A. Snipes and the ITER H-mode Threshold Database Working Group, "An Analysis of the H-mode Threshold in ITER,"
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
        - A scaling with the total ion mass is used in this model. Martin 08 shows that P_LH scales with 1/m_i. It is stated;
        "When this mass dependence is applied to the deuterium-tritium discharges for ITER, the above predicted values of P_LH can be
        reduced by ~ 20%". We thus apply a (2/m_i) addition so that for a 50/50 D-T mixture (M_i = 2.5 amu), the predicted values is 20% lower.


    References
    ----------
        - Y. R. Martin, T. Takizuka, and the I. C. H-mode. T. D. Group, “Power requirement for accessing the H-mode in ITER,”
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
        - A scaling with the total ion mass is used in this model. Martin 08 shows that P_LH scales with 1/m_i. It is stated;
        "When this mass dependence is applied to the deuterium-tritium discharges for ITER, the above predicted values of P_LH can be
        reduced by ~ 20%". We thus apply a (2/m_i) addition so that for a 50/50 D-T mixture (M_i = 2.5 amu), the predicted values is 20% lower.
    References
    ----------
        - Y. R. Martin, T. Takizuka, and the I. C. H-mode. T. D. Group, “Power requirement for accessing the H-mode in ITER,”
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
        - A scaling with the total ion mass is used in this model. Martin 08 shows that P_LH scales with 1/m_i. It is stated;
        "When this mass dependence is applied to the deuterium-tritium discharges for ITER, the above predicted values of P_LH can be
        reduced by ~ 20%". We thus apply a (2/m_i) addition so that for a 50/50 D-T mixture (M_i = 2.5 amu), the predicted values is 20% lower.

    References
    ----------
        - Y. R. Martin, T. Takizuka, and the I. C. H-mode. T. D. Group, “Power requirement for accessing the H-mode in ITER,”
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
        - A scaling with the total ion mass is used in this model. Snipes cites that P_LH scales with 1/m_i. It is stated;
        "This results in a 20% reduction in the threshold power for a 50/50 D-T mixture compared with the pure deuterium
        results above". We thus apply a (2/m_i) addition so that for a 50/50 D-T mixture (M_i = 2.5 amu), the predicted values is 20% lower.

    References
    ----------
        - J. A. Snipes and the I. H-mode. T. Group, “Latest results on the H-mode threshold using the international H-mode threshold database,”
        Plasma Physics and Controlled Fusion, vol. 42, no. 5A, pp. A299-A308, May 2000,
        doi: https://doi.org/10.1088/0741-3335/42/5a/336.
    """
    return (
        1.42
        * dnla20**0.58
        * b_plasma_toroidal_on_axis**0.82
        * rmajor
        * rminor**0.81
        * (2.0 / m_ions_total_amu)
    )


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
        - A scaling with the total ion mass is used in this model. Snipes cites that P_LH scales with 1/m_i. It is stated;
        "This results in a 20% reduction in the threshold power for a 50/50 D-T mixture compared with the pure deuterium
        results above". We thus apply a (2/m_i) addition so that for a 50/50 D-T mixture (M_i = 2.5 amu), the predicted values is 20% lower.

    References
    ----------
        - J. A. Snipes and the I. H-mode. T. Group, “Latest results on the H-mode threshold using the international H-mode threshold database,”
        Plasma Physics and Controlled Fusion, vol. 42, no. 5A, pp. A299-A308, May 2000,
        doi: https://doi.org/10.1088/0741-3335/42/5a/336.

    """
    return (
        1.547
        * dnla20**0.615
        * b_plasma_toroidal_on_axis**0.851
        * rmajor**1.089
        * rminor**0.876
        * (2.0 / m_ions_total_amu)
    )


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
        - A scaling with the total ion mass is used in this model. Snipes cites that P_LH scales with 1/m_i. It is stated;
        "This results in a 20% reduction in the threshold power for a 50/50 D-T mixture compared with the pure deuterium
        results above". We thus apply a (2/m_i) addition so that for a 50/50 D-T mixture (M_i = 2.5 amu), the predicted values is 20% lower.

    References
    ----------
        - J. A. Snipes and the I. H-mode. T. Group, “Latest results on the H-mode threshold using the international H-mode threshold database,”
        Plasma Physics and Controlled Fusion, vol. 42, no. 5A, pp. A299-A308, May 2000,
        doi: https://doi.org/10.1088/0741-3335/42/5a/336.

    """
    return (
        1.293
        * dnla20**0.545
        * b_plasma_toroidal_on_axis**0.789
        * rmajor**0.911
        * rminor**0.744
        * (2.0 / m_ions_total_amu)
    )


def calculate_snipes2000_closed_divertor_nominal(
    dnla20: float,
    b_plasma_toroidal_on_axis: float,
    rmajor: float,
    m_ions_total_amu: float,
) -> float:
    """Calculate the nominal Snipes 2000 Closed Divertor L-H transition power threshold with CD factor.

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
        - A scaling with the total ion mass is used in this model. Snipes cites that P_LH scales with 1/m_i. It is stated;
        "This results in a 20% reduction in the threshold power for a 50/50 D-T mixture compared with the pure deuterium
        results above". We thus apply a (2/m_i) addition so that for a 50/50 D-T mixture (M_i = 2.5 amu), the predicted values is 20% lower.

    References
    ----------
        - J. A. Snipes and the I. H-mode. T. Group, “Latest results on the H-mode threshold using the international H-mode threshold database,”
        Plasma Physics and Controlled Fusion, vol. 42, no. 5A, pp. A299-A308, May 2000,
        doi: https://doi.org/10.1088/0741-3335/42/5a/336.

    """
    return (
        0.8
        * dnla20**0.5
        * b_plasma_toroidal_on_axis**0.53
        * rmajor**1.51
        * (2.0 / m_ions_total_amu)
    )


def calculate_snipes2000_closed_divertor_upper(
    dnla20: float,
    b_plasma_toroidal_on_axis: float,
    rmajor: float,
    m_ions_total_amu: float,
) -> float:
    """Calculate the upper Snipes 2000 Closed Divertor L-H transition power threshold with CD factor.

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
        - A scaling with the total ion mass is used in this model. Snipes cites that P_LH scales with 1/m_i. It is stated;
        "This results in a 20% reduction in the threshold power for a 50/50 D-T mixture compared with the pure deuterium
        results above". We thus apply a (2/m_i) addition so that for a 50/50 D-T mixture (M_i = 2.5 amu), the predicted values is 20% lower.

    References
    ----------
        - J. A. Snipes and the I. H-mode. T. Group, “Latest results on the H-mode threshold using the international H-mode threshold database,”
        Plasma Physics and Controlled Fusion, vol. 42, no. 5A, pp. A299-A308, May 2000,
        doi: https://doi.org/10.1088/0741-3335/42/5a/336.

    """
    return (
        0.867
        * dnla20**0.561
        * b_plasma_toroidal_on_axis**0.588
        * rmajor**1.587
        * (2.0 / m_ions_total_amu)
    )


def calculate_snipes2000_closed_divertor_lower(
    dnla20: float,
    b_plasma_toroidal_on_axis: float,
    rmajor: float,
    m_ions_total_amu: float,
) -> float:
    """Calculate the lower Snipes 2000 Closed Divertor L-H transition power threshold with CD factor.

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
        - A scaling with the total ion mass is used in this model. Snipes cites that P_LH scales with 1/m_i. It is stated;
        "This results in a 20% reduction in the threshold power for a 50/50 D-T mixture compared with the pure deuterium
        results above". We thus apply a (2/m_i) addition so that for a 50/50 D-T mixture (M_i = 2.5 amu), the predicted values is 20% lower.

    References
    ----------
        - J. A. Snipes and the I. H-mode. T. Group, “Latest results on the H-mode threshold using the international H-mode threshold database,”
        Plasma Physics and Controlled Fusion, vol. 42, no. 5A, pp. A299-A308, May 2000,
        doi: https://doi.org/10.1088/0741-3335/42/5a/336.

    """
    return (
        0.733
        * dnla20**0.439
        * b_plasma_toroidal_on_axis**0.472
        * rmajor**1.433
        * (2.0 / m_ions_total_amu)
    )


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

    Notes
    -----

    References
    ----------
        - A. E. Hubbard et al., “Threshold conditions for transitions to I-mode and H-mode with unfavourable ion grad B drift direction,”
        Nuclear Fusion, vol. 52, no. 11, pp. 114009-114009, Oct. 2012,
        doi: https://doi.org/10.1088/0029-5515/52/11/114009.

    """
    return 2.11 * (plasma_current / 1e6) ** 0.94 * dnla20**0.65


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

    Notes
    -----

    References
    ----------
        - A. E. Hubbard et al., “Threshold conditions for transitions to I-mode and H-mode with unfavourable ion grad B drift direction,”
        Nuclear Fusion, vol. 52, no. 11, pp. 114009-114009, Oct. 2012,
        doi: https://doi.org/10.1088/0029-5515/52/11/114009.

    """
    return 2.11 * (plasma_current / 1e6) ** 1.18 * dnla20**0.83


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

    Notes
    -----

    References
    ----------
        - A. E. Hubbard et al., “Threshold conditions for transitions to I-mode and H-mode with unfavourable ion grad B drift direction,”
        Nuclear Fusion, vol. 52, no. 11, pp. 114009-114009, Oct. 2012,
        doi: https://doi.org/10.1088/0029-5515/52/11/114009.

    """
    return 2.11 * (plasma_current / 1e6) ** 0.7 * dnla20**0.47


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
        - A. E. Hubbard et al., “Physics and performance of the I-mode regime over an expanded operating space on Alcator C-Mod,”
        Nuclear Fusion, vol. 57, no. 12, p. 126039, Oct. 2017,
        doi: https://doi.org/10.1088/1741-4326/aa8570.

    """
    return 0.162 * dnla20 * a_plasma_surface * b_plasma_toroidal_on_axis**0.26


def calculate_martin08_aspect_nominal(
    dnla20: float,
    b_plasma_toroidal_on_axis: float,
    a_plasma_surface: float,
    m_ions_total_amu: float,
    aspect: float,
) -> float:
    """Calculate the nominal Martin L-H transition power threshold with aspect ratio correction from T Takizuka.

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
        - Thus will return an aspect ratio correction of the aspect ratio is less than or equal to 2.7.
        if not the usual Martin 2008 scaling will be returned.

        - A scaling with the total ion mass is used in this model. Martin 08 shows that P_LH scales with 1/m_i. It is stated;
        "When this mass dependence is applied to the deuterium-tritium discharges for ITER, the above predicted values of P_LH can be
        reduced by ~ 20%". We thus apply a (2/m_i) addition so that for a 50/50 D-T mixture (M_i = 2.5 amu), the predicted values is 20% lower.

    References
    ----------
        - Y. R. Martin, T. Takizuka, and the I. C. H-mode. T. D. Group, “Power requirement for accessing the H-mode in ITER,”
        Journal of Physics: Conference Series, vol. 123, p. 012033, Jul. 2008,
        doi: https://doi.org/10.1088/1742-6596/123/1/012033.

        - T. Takizuka et.al, “Roles of aspect ratio, absolute B and effective Z of the H-mode power threshold in tokamaks of the ITPA database,”
        Plasma Physics and Controlled Fusion, vol. 46, no. 5A, pp. A227-A233, Apr. 2004,
        doi: https://doi.org/10.1088/0741-3335/46/5a/024.

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


def calculate_martin08_aspect_upper(
    dnla20: float,
    b_plasma_toroidal_on_axis: float,
    a_plasma_surface: float,
    m_ions_total_amu: float,
    aspect: float,
) -> float:
    """Calculate the upper Martin L-H transition power threshold with aspect ratio correction from T Takizuka.

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
        - Thus will return an aspect ratio correction of the aspect ratio is less than or equal to 2.7.
        if not the usual Martin 2008 scaling will be returned.

        - A scaling with the total ion mass is used in this model. Martin 08 shows that P_LH scales with 1/m_i. It is stated;
        "When this mass dependence is applied to the deuterium-tritium discharges for ITER, the above predicted values of P_LH can be
        reduced by ~ 20%". We thus apply a (2/m_i) addition so that for a 50/50 D-T mixture (M_i = 2.5 amu), the predicted values is 20% lower.

    References
    ----------
        - Y. R. Martin, T. Takizuka, and the I. C. H-mode. T. D. Group, “Power requirement for accessing the H-mode in ITER,”
        Journal of Physics: Conference Series, vol. 123, p. 012033, Jul. 2008,
        doi: https://doi.org/10.1088/1742-6596/123/1/012033.

        - T. Takizuka et.al, “Roles of aspect ratio, absolute B and effective Z of the H-mode power threshold in tokamaks of the ITPA database,”
        Plasma Physics and Controlled Fusion, vol. 46, no. 5A, pp. A227-A233, Apr. 2004,
        doi: https://doi.org/10.1088/0741-3335/46/5a/024.

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


def calculate_martin08_aspect_lower(
    dnla20: float,
    b_plasma_toroidal_on_axis: float,
    a_plasma_surface: float,
    m_ions_total_amu: float,
    aspect: float,
) -> float:
    """Calculate the lower Martin L-H transition power threshold with aspect ratio correction from T Takizuka.

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
        - Thus will return an aspect ratio correction of the aspect ratio is less than or equal to 2.7.
        if not the usual Martin 2008 scaling will be returned.

        - A scaling with the total ion mass is used in this model. Martin 08 shows that P_LH scales with 1/m_i. It is stated;
        "When this mass dependence is applied to the deuterium-tritium discharges for ITER, the above predicted values of P_LH can be
        reduced by ~ 20%". We thus apply a (2/m_i) addition so that for a 50/50 D-T mixture (M_i = 2.5 amu), the predicted values is 20% lower.

    References
    ----------
        - Y. R. Martin, T. Takizuka, and the I. C. H-mode. T. D. Group, “Power requirement for accessing the H-mode in ITER,”
        Journal of Physics: Conference Series, vol. 123, p. 012033, Jul. 2008,
        doi: https://doi.org/10.1088/1742-6596/123/1/012033.

        - T. Takizuka et.al, “Roles of aspect ratio, absolute B and effective Z of the H-mode power threshold in tokamaks of the ITPA database,”
        Plasma Physics and Controlled Fusion, vol. 46, no. 5A, pp. A227-A233, Apr. 2004,
        doi: https://doi.org/10.1088/0741-3335/46/5a/024.

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
