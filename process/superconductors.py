import logging

import numpy as np
from scipy import optimize

from process.exceptions import ProcessValueError
from process.fortran import error_handling as eh
from process.fortran import rebco_variables

logger = logging.getLogger(__name__)


def jcrit_rebco(temperature, b):
    """Critical current density for "REBCO" 2nd generation HTS superconductor
    temperature : input real : superconductor temperature (K)
    b : input real : Magnetic field at superconductor (T)
    jcrit : output real : Critical current density in superconductor (A/m2)

    Will return a negative number if the temperature is greater than Tc0, the
    zero-field critical temperature.
    """
    tc0 = 90.0  # (K)
    birr0 = 132.5  # (T)
    a = 1.82962e8  # scaling constant
    # exponents
    p = 0.5875
    q = 1.7
    alpha = 1.54121
    beta = 1.96679
    oneoveralpha = 1 / alpha

    validity = True

    if (temperature < 4.2) or (temperature > 72.0):
        validity = False
    if temperature < 65:
        if (b < 0.0) or (b > 15.0):
            validity = False
    else:
        if (b < 0.0) or (b > 11.5):
            validity = False

    if not validity:
        logger.warning(
            f"""jcrit_rebco: input out of range
            temperature: {temperature}
            Field: {b}
            """
        )

    if temperature < tc0:
        # Normal case
        birr = birr0 * (1 - temperature / tc0) ** alpha
    else:
        # If temp is greater than critical temp, ensure result is real but negative.
        birr = birr0 * (1 - temperature / tc0)

    if b < birr:
        # Normal case
        factor = (b / birr) ** p * (1 - b / birr) ** q
        jcrit = (a / b) * (birr**beta) * factor
    else:
        # Field is too high
        # Ensure result is real but negative, and varies with temperature.
        # tcb = critical temperature at field b
        tcb = tc0 * (1 - (b / birr0) ** oneoveralpha)
        jcrit = -(temperature - tcb)

    return jcrit, validity


def current_sharing_rebco(bfield, j):
    """Current sharing temperature for "REBCO" 2nd generation HTS superconductor
    b : input real : Magnetic field at superconductor (T)
    j : input real : Current density in superconductor (A/m2)
    current_sharing_t : output real : Current sharing temperature (K)
    """

    def deltaj_rebco(temperature):
        jcritical, _ = jcrit_rebco(temperature, bfield)
        return jcritical - j

    # No additional arguments are required for deltaj_rebco since it only has one argument.

    estimate = 10.0
    another_estimate = 20.0
    current_sharing_t, root_result = optimize.newton(
        deltaj_rebco,
        estimate,
        tol=1e-6,
        rtol=1e-6,
        maxiter=50,
        x1=another_estimate,
        full_output=True,
        disp=True,
    )

    return current_sharing_t


def itersc(
    temp_conductor: float,
    b_conductor: float,
    strain: float,
    b_c20max: float,
    temp_c0max: float,
) -> tuple[float, float, float]:
    """
    Calculate the critical current density, critical field, and critical temperature
    for an ITER Nb₃Sn superconductor using the ITER Nb₃Sn critical surface model.

    :param temp_conductor: Superconductor temperature (K).
    :type temp_conductor: float
    :param b_conductor: Magnetic field at the conductor (T).
    :type b_conductor: float
    :param strain: Strain in the superconductor.
    :type strain: float
    :param b_c20max: Upper critical field (T) for the superconductor at zero temperature and strain.
    :type b_c20max: float
    :param temp_c0max: Critical temperature (K) at zero field and strain.
    :type temp_c0max: float
    :return: A tuple containing:
        - j_critical: Critical current density in the superconductor (A/m²).
        - b_critical: Critical field (T).
        - temp_critical: Critical temperature (K).
    :rtype: tuple[float, float, float]

    :notes:
        - This routine uses the ITER Nb₃Sn critical surface model.
        - The model assumes a strand size of 0.82 mm diameter.

    :references:
        - ITER DDD 11-7: Magnets - conductors (2NBKXY) (2009),
          https://user.iter.org/?uid=2NBKXY&action=get_document
    """
    # Scaling constant C [AT/mm²]
    csc = 19922.0
    # Low field exponent
    p = 0.63
    # High field exponent
    q = 2.1
    # Strain fitting constant C_{a1}
    ca1 = 44.48
    # Strain fitting constant C_{a2}
    ca2 = 0.0
    # Residual strain component epsilon_{0,a}
    epsilon_0a = 0.00256

    # ITER strand diameter (mm)
    diter = 0.82

    # ITER strand copper fraction
    f_a_strand_copper = 0.5

    j_scaling, b_critical, temp_critical = bottura_scaling(
        csc=csc,
        p=p,
        q=q,
        c_a1=ca1,
        c_a2=ca2,
        epsilon_0a=epsilon_0a,
        temp_conductor=temp_conductor,
        b_conductor=b_conductor,
        epsilon=strain,
        b_c20max=b_c20max,
        temp_c0max=temp_c0max,
    )

    # Critical current density in superconductor (A/m²)
    # ITER parameters are for the current in a single strand,
    # not per unit area, so scalefac converts to current density.
    # Convert from mm² to m².
    scalefac = np.pi * (0.5 * diter) ** 2 * (1.0 - f_a_strand_copper) / 1.0e6

    j_critical = j_scaling / scalefac

    return j_critical, b_critical, temp_critical


def jcrit_nbti(
    temp_conductor: float,
    b_conductor: float,
    c0: float,
    b_c20max: float,
    temp_c0max: float,
) -> tuple[float, float]:
    """
        Calculate the critical current density and critical temperature for a NbTi superconductor strand using the old empirical
        Lubell scaling law.

        :param temp_conductor: Superconductor temperature (K).
        :type temp_conductor: float
        :param b_conductor: Magnetic field at the conductor (T).
        :type b_conductor: float
        :param c0: Scaling constant (A/m²).
        :type c0: float
        :param b_c20max: Upper critical field (T) for the superconductor at zero temperature and strain.
        :type b_c20max: float
        :param temp_c0max: Critical temperature (K) at zero field and strain.
        :type temp_c0max: float
        :return: A tuple containing:
            - jcrit: Critical current density in the superconductor (A/m²).
            - tcrit: Critical temperature (K).
        :rtype: tuple[float, float]

        :notes:
            - If the magnetic field exceeds the upper critical field (bmax > bc20max),
              the critical temperature is adjusted to ensure a real (negative) value.
            - If the temperature exceeds the critical temperature, the critical surface
              is considered exceeded, and the reduced temperature (tbar) becomes negative.

            - This model uses an antiquated scaling law for NbTi, which is not used in the superconductor field for nearly 30 years.
            - It is simplistic and linear in J_c (B) (for a fixed temperature), lacking accuracy in both the high and low field regions.

        :references:
            - M. Lubell, “Empirical scaling formulas for critical current and critical field for commercial NbTi,”
            IEEE Transactions on Magnetics, vol. 19, no. 3, pp. 754-757, May 1983,
            doi: https://doi.org/10.1109/tmag.1983.1062311.
    ‌
    """
    bratio = b_conductor / b_c20max

    # Critical temperature (K) at field
    temp_critical = (
        (temp_c0max * (1.0 - bratio) ** 0.59)
        if bratio < 1
        else (temp_c0max * (1.0 - bratio))
    )

    # Reduced temperature
    tbar = 1.0 - temp_conductor / temp_critical

    # Critical current density (A/m²)
    j_critical = c0 * (1.0 - bratio) * tbar

    return j_critical, temp_critical


def bi2212(b_conductor, jstrand, temp_conductor, f_strain):
    """
        Fitted parameterization to Bi-2212 superconductor properties.

        This function calculates the critical current density and the temperature margin
        for Bi-2212 superconductor in the TF coils using a fit by M. Kovari to measurements
        described in the reference, specifically from the points shown in Figure 6.

        Bi-2212 (Bi₂Sr₂CaCu₂O₈₋ₓ) is a first-generation high-temperature superconductor.
        It needs to be operated below about 10K but remains superconducting at much higher
        fields at that temperature than Nb₃Sn, etc.

        :param b_conductor: Magnetic field at conductor (T).
        :type b_conductor: float
        :param jstrand: Current density in strand (A/m²).
        :type jstrand: float
        :param temp_conductor: Superconductor temperature (K).
        :type temp_conductor: float
        :param f_strain: Adjustment factor (≤ 1) to account for strain, radiation damage,
                     fatigue, or AC losses.
        :type f_strain: float
        :raises ProcessValueError: If the input parameters are outside the range of validity.
        :return: A tuple containing:
            - j_critical: Critical current density in strand (A/m²).
            - temp_margin: Temperature margin (K).
        :rtype: tuple[float, float]

        :notes:
            -The model's range of validity is:
                T < 20K
                Adjusted field b < 104 T
                B > 6 T

        :reference:
            - D. C. Larbalestier, J. Jiang, U. P. Trociewitz, F. Kametani, and E. E. Hellstrom,
            “A transformative superconducting magnet technology for fields well above 30 T using isotropic round wire multifilament Bi2Sr2CaCu2O8-x conductor,”
            May 06, 2013. https://www.researchgate.net/publication/236627864_A_transformative_superconducting_magnet_technology_for_fields_well_above_30_T_using_isotropic_round_wire_multifilament_Bi2Sr2CaCu2O8-x_conductor
    ‌
    """

    b = b_conductor / np.exp(-0.168 * (temp_conductor - 4.2))

    #  Engineering (i.e. strand) critical current density (A/m2)

    j_critical = f_strain * (1.175e9 * np.exp(-0.02115 * b) - 1.288e8)

    #  Temperature margin (K)
    #  Simple inversion of above calculation, using actual current density
    #  in strand instead of jcrit

    temp_margin = (
        1.0
        / 0.168
        * np.log(
            np.log(1.175e9 / (jstrand / f_strain + 1.288e8)) / (0.02115 * b_conductor)
        )
        + 4.2
        - temp_conductor
    )

    #  Check if ranges of validity have been violated

    if (temp_conductor > 20.0) or (b_conductor < 6.0) or (b > 104.0):
        raise ProcessValueError(
            "Fit extrapolated outside of range of validity",
            temp_conductor=temp_conductor,
            b_conductor=b_conductor,
            b=b,
        )

    return j_critical, temp_margin


def gl_nbti(
    temp_conductor: float,
    b_conductor: float,
    strain: float,
    b_c20max: float,
    t_c0: float,
) -> tuple[float, float, float]:
    """
    Calculate the critical current density, critical field, and critical temperature
    of an ITER Nb-Ti strand based on the Ginzburg-Landau theory of superconductivity.

    :param temp_conductor: Superconductor temperature [K].
    :type temp_conductor: float
    :param b_conductor: Magnetic field at the conductor [T].
    :type b_conductor: float
    :param strain: Intrinsic strain in the superconductor [%].
    :type strain: float
    :param b_c20max: Strain-dependent upper critical field at zero temperature [T].
    :type b_c20max: float
    :param t_c0: Strain-dependent critical temperature at zero strain [K].
    :type t_c0: float

    :return: A tuple containing:
        - j_critical : Critical current density in the superconductor [A/m²].
        - b_critical : Critical magnetic field [T].
        - t_critical : Critical temperature [K].
    :rtype: tuple[float, float, float]

    :notes:


    :references:

    """

    a_0 = 1102e6
    p = 0.49
    q = 0.56
    n = 1.83
    v = 1.42
    u = 0.0
    w = 2.2

    # Strain fitted parameters to a single filament
    c2 = -0.0025
    c3 = -0.0003
    c4 = -0.0001
    epsilon_m = -0.002e-2

    # ==========================================================

    epsilon_i = strain - epsilon_m

    strain_func = (
        1 + c2 * (epsilon_i) ** 2 + c3 * (epsilon_i) ** 3 + c4 * (epsilon_i) ** 4
    )

    t_e = t_c0 * strain_func ** (1 / w)

    t_reduced = temp_conductor / t_e

    a_e = a_0 * strain_func ** (u / w)

    # Critical Field
    b_critical = b_c20max * (1 - t_reduced**v) * strain_func

    b_reduced = b_conductor / b_critical

    # Critical temperature (K)
    t_critical = t_e

    # Critical current density (A/m2)
    if b_reduced <= 1.0:
        j_critical = (
            a_e
            * (t_e * (1 - t_reduced**2)) ** 2
            * b_critical ** (n - 3)
            * b_reduced ** (p - 1)
            * (1 - b_reduced) ** q
        )
    else:  # Fudge to yield negative single valued function of Jc for fields above Bc2
        j_critical = (
            a_e
            * (t_e * (1 - t_reduced**2)) ** 2
            * b_critical ** (n - 3)
            * b_reduced ** (p - 1)
            * (1 - b_reduced) ** 1.0
        )

    return j_critical, b_critical, t_critical


def gl_rebco(thelium, bmax, strain, bc20max, t_c0):
    """Author: S B L Chislett-McDonald Durham University
    Category: subroutine

    Critical current density of a SuperPower REBCO tape based on measurements by P. Branch
    at Durham University
    https://git.ccfe.ac.uk/process/process/uploads/e98c6ea13da782cdc6fe16daea92078a/20200707_Branch-Osamura-Hampshire_-_accepted_SuST.pdf
    and fit to state-of-the-art measurements at 4.2 K published in SuST
    http://dx.doi.org/10.1088/0953-2048/24/3/035001

    \\begin{equation}
    J_{c,TS}(B,T,\\epsilon_{I}) = A(\\epsilon_{I}) \\left[T_{c}(\\epsilon_{I})*(1-t^2)\\right]^2\\left
    [B_{c2}(\\epsilon_I)*(1-t)^s\\right]^{n-3}b^{p-1}(1-b)^q~.
    \\end{equation}

    - \\( \\thelium \\) -- Coolant/SC temperature [K]
    - \\( \\bmax \\) -- Magnetic field at conductor [T]
    - \\( \\epsilon_{I} \\) -- Intrinsic strain in superconductor [\\%]
    - \\( \\B_{c2}(\\epsilon_I) \\) -- Strain dependent upper critical field [T]
    - \\( \\b \\) -- Reduced field = bmax / \\B_{c2}(\\epsilon_I)*(1-t^\\nu) [unitless]
    - \\( \\T_{c}(\\epsilon_{I}) \\) -- Strain dependent critical temperature (K)
    - \\( \\t \\) -- Reduced temperature = thelium / \\T_{c}(\\epsilon_{I}) [unitless]
    - \\( \\A(epsilon_{I}) \\) -- Strain dependent Prefactor [A / ( m\\(^2\\) K\\(^-2) T\\(^n-3))]
    - \\( \\J_{c,TS} \\) --  Critical current density in superconductor [A / m\\(^-2\\)]
    - \\( \\epsilon_{m} \\) -- Strain at which peak in J_c occurs [\\%]
    """
    # critical current density prefactor
    a_0 = 2.95e2
    # flux pinning field scaling parameters
    p = 0.32
    q = 2.50
    n = 3.33
    # temperatute scaling parameter
    s = 5.27
    # strain scaling parameters
    c2 = -0.0191
    c3 = 0.0039
    c4 = 0.00103
    em = 0.058
    # strain conversion parameters
    u = 0.0
    w = 2.2

    epsilon_i = strain - em

    strain_func = (
        1 + c2 * (epsilon_i) ** 2 + c3 * (epsilon_i) ** 3 + c4 * (epsilon_i) ** 4
    )

    t_e = t_c0 * strain_func ** (1 / w)

    t_reduced = thelium / t_e

    a_e = a_0 * strain_func ** (u / w)

    #  Critical Field
    bcrit = bc20max * (1 - t_reduced) ** s * strain_func

    b_reduced = bmax / bcrit

    #  Critical temperature (K)
    tcrit = t_e

    #  Critical current density (A/m2)
    jcrit = (
        a_e
        * (t_e * (1 - t_reduced**2)) ** 2
        * bcrit ** (n - 3)
        * b_reduced ** (p - 1)
        * (1 - b_reduced) ** q
    )

    return jcrit, bcrit, tcrit


def hijc_rebco(thelium, bmax, _strain, bc20max, t_c0):
    """Implementation of High Current Density REBCO tape
    author: R Chapman, UKAEA
    thelium : input real : SC temperature (K)
    bmax : input real : Magnetic field at conductor (T)
    strain : input real : Strain in superconductor
    bc20max : input real : Upper critical field (T) for superconductor
    at zero temperature and strain
    t_c0 : input real : Critical temperature (K) at zero field and strain
    jcrit : output real : Critical current density in superconductor (A/m2)
    bcrit : output real : Critical field (T)
    tcrit : output real : Critical temperature (K)

    Returns the critical current of a REBCO tape based on a critical surface
    (field, temperature) parameterization. Based in part on the parameterization
    described in: M. J. Wolf, N. Bagrets, W. H. Fietz, C. Lange and K. Weiss,
    "Critical Current Densities of 482 A/mm2 in HTS CrossConductors at 4.2 K and 12 T,"
    in IEEE Transactions on Applied Superconductivity, vol. 28, no. 4, pp. 1-4,
    June 2018, Art no. 4802404, doi: 10.1109/TASC.2018.2815767.

    And on the experimental
    data presented here: "2G HTS Wire Development at SuperPower", Drew W. Hazelton,
    February 16, 2017 https://indico.cern.ch/event/588810/contributions/2473740/
    The high Ic parameterization is a result of modifications based on Ic values
    observed in: "Conceptual design of HTS magnets for fusion nuclear science facility",
    Yuhu Zhai, Danko van der Laan, Patrick Connolly, Charles Kessel, 2021,
    https://doi.org/10.1016/j.fusengdes.2021.112611

    The parameter A is transformed into a function A(T) based on a Newton polynomial fit
    considering A(4.2 K) = 2.2e8, A(20 K) = 2.3e8 and A(65 K) = 3.5e8. These values were
    selected manually. A good fit to the pubished data can be seen in the 4-10 T range
    but the fit deviates at very low or very high field.

    C. Ashe 2/11/23
    A scaling factor of 0.4 was originally applied to jcrit to accomodate the models
    original purpose of being used with CORC cables, which PROCESS did not model as it
    assumed a stacked tape block design.
    """

    a = 1.4
    b = 2.005
    # critical current density prefactor
    a_0 = 2.2e8
    # flux pinning field scaling parameters
    p = 0.39
    q = 0.9
    # strain conversion parameters
    u = 33450.0
    v = -176577.0

    # Critical Field (T)
    # B_crit(T) calculated using temperature and critical temperature
    bcrit = bc20max * (1.0 - thelium / t_c0) ** a

    # Critical temperature (K)
    # scaled to match behaviour in GL_REBCO routine,
    # ONLY TO BE USED until a better suggestion is received
    tcrit = 0.999965 * t_c0

    # finding A(T); constants based on a Newton polynomial fit to pubished data
    a_t = a_0 + (u * thelium**2) + (v * thelium)

    # Critical current density (A/m2)
    # In the original formula bcrit must be > bmax to prevent NaNs.
    # However, negative jcrit is permissible (I think).
    # So when bcrit < bmax, I reverse the sign of the bracket,
    # giving a negative but real value of jcrit.

    if bcrit > bmax:
        jcrit = (a_t / bmax) * bcrit**b * (bmax / bcrit) ** p * (1 - bmax / bcrit) ** q
    else:
        jcrit = (a_t / bmax) * bcrit**b * (bmax / bcrit) ** p * (bmax / bcrit - 1) ** q

    # print("thelium = ", thelium, "   bcrit = ", bcrit, "   bmax = ", bmax, "   1 - bmax / bcrit = ", 1 - bmax / bcrit)

    # Jc times HTS area: default area is width 4mm times HTS layer thickness 1 um,
    # divided by the tape area to provide engineering Jc per tape,!
    # A scaling factor of 0.4 used to be applied below to assume the difference
    # between tape stacks and CORC cable layouts.
    jcrit = (
        jcrit
        * (rebco_variables.tape_width * rebco_variables.rebco_thickness)
        / (rebco_variables.tape_width * rebco_variables.tape_thickness)
    )

    return jcrit, bcrit, tcrit


def western_superconducting_nb3sn(
    temp_conductor: float,
    b_conductor: float,
    strain: float,
    b_c20max: float,
    temp_c0max: float,
) -> tuple[float, float, float]:
    """
    Calculate the critical current density, critical field, and critical temperature
    for a WST Nb₃Sn superconductor using the ITER Nb₃Sn critical surface model.

    :param temp_conductor: Superconductor temperature (K).
    :type temp_conductor: float
    :param b_conductor: Magnetic field at the superconductor (T).
    :type b_conductor: float
    :param strain: Strain in the superconductor.
    :type strain: float
    :param b_c20max: Upper critical field (T) for the superconductor at zero temperature and strain.
    :type b_c20max: float
    :param temp_c0max: Critical temperature (K) at zero field and strain.
    :type temp_c0max: float
    :return: A tuple containing:
        - j_critical: Critical current density in the superconductor (A/m²).
        - b_crititical: Critical field (T).
        - temp_critical: Critical temperature (K).
    :rtype: tuple[float, float, float]

    :notes:
        - This routine uses the WST Nb3Sn critical surface model.
        - The scaling constants and parameters are based on the reference below.
        - This assumes a strand size of 1.5mm.
        - Compared to the EUTF4 (OST) ( European qualification samples for TF conductor),
        the performance of the WST at low strain is superior by about 10%.

    :references:
        - V. Corato, “EUROFUSION WPMAG-REP(16) 16565 Common operating values for DEMO magnets design for 2016 REPORT.”
        Accessed: May 12, 2025. [Online].
        Available: https://scipub.euro-fusion.org/wp-content/uploads/eurofusion/WPMAGREP16_16565_submitted.pdf

        - “Introduction of WST,” 2015. Accessed: May 12, 2025. [Online].
        Available: https://indico.cern.ch/event/340703/contributions/802232/attachments/668814/919331/WST_INTRO_2015-3_for_FCC_WEEK.pdf

    """
    # Scaling constant C [AT/mm²]
    csc = 83075.0
    # Low field exponent p
    p = 0.593
    # High field exponent q
    q = 2.156
    # Strain fitting constant C_{a1}
    c_a1 = 50.06
    # Strain fitting constant C_{a2}
    c_a2 = 0.0
    # Residual strain component epsilon_{0,a}
    epsilon_0a = 0.00312

    j_scaling, b_critical, t_critical = bottura_scaling(
        csc=csc,
        p=p,
        q=q,
        c_a1=c_a1,
        c_a2=c_a2,
        epsilon_0a=epsilon_0a,
        temp_conductor=temp_conductor,
        b_conductor=b_conductor,
        epsilon=strain,
        b_c20max=b_c20max,
        temp_c0max=temp_c0max,
    )

    # Scale from mm² to m²
    scalefac = 1.0e6
    j_critical = j_scaling * scalefac
    return j_critical, b_critical, t_critical


def bottura_scaling(
    csc: float,
    p: float,
    q: float,
    c_a1: float,
    c_a2: float,
    epsilon_0a: float,
    temp_conductor: float,
    b_conductor: float,
    epsilon: float,
    b_c20max: float,
    temp_c0max: float,
) -> tuple[float, float, float]:
    """
    Implements the scaling from:
    Jc(B,T,epsilon) Parameterization for the ITER Nb₃Sn Production,

    :param csc: Scaling constant C [AT/mm²].
    :type csc: float
    :param p: Low field exponent of the pinning force
    :type p: float
    :param q: High field exponent of the pinning force
    :type q: float
    :param c_a1: Strain fitting constant C_{a1}.
    :type c_a1: float
    :param c_a2: Strain fitting constant C_{a2}.
    :type c_a2: float
    :param epsilon_0a: Residual strain component
    :type epsilon_0a: float
    :param temp_conductor: Superconductor temperature (K).
    :type temp_conductor: float
    :param b_conductor: Magnetic field at conductor (T).
    :type b_conductor: float
    :param epsilon: Strain in superconductor.
    :type epsilon: float
    :param b_c20max: Upper critical field (T) for superconductor at zero temperature and strain.
    :type b_c20max: float
    :param temp_c0max: Critical temperature (K) at zero field and strain.
    :type temp_c0max: float
    :return: A tuple containing:
        - jscaling: Critical current density scaling factor (A/m²).
        - bcrit: Critical field (T).
        - tcrit: Critical temperature (K).
    :rtype: tuple[float, float, float]

    :notes:
        - This is a generic scaling proposed for the characterization and production of
        ITER Nb₃Sn strands. This is also known as the "ITER-2008 parametrization."

        - Parameter ranges are strain (1.5% to 0.4%), temperature (2.35 to 16 K), and field (0.5 to 19 T).
        The ITER-2008 parameterization achieves an average accuracy error of 3.8 Amps, with the best at 1.5 Amps and the worst at 7.5 Amps.

        - The strain function is suitable only in the moderate strain region, down to 0.8%.

    :references:
        - L. Bottura and B. Bordini, “$J_{C}(B,T,\varepsilon)$ Parameterization for the ITER ${\rm Nb}_{3}{\rm Sn}$ Production,”
        IEEE Transactions on Applied Superconductivity, vol. 19, no. 3, pp. 1521-1524, Jun. 2009,
        doi: https://doi.org/10.1109/tasc.2009.2018278.
    """

    epsilon_sh = (c_a2 * epsilon_0a) / (np.sqrt(c_a1**2 - c_a2**2))

    # Strain function
    # 0.83 < s < 1.0, for -0.005 < strain < 0.005
    strfun = np.sqrt(epsilon_sh**2 + epsilon_0a**2) - np.sqrt(
        (epsilon - epsilon_sh) ** 2 + epsilon_0a**2
    )
    strfun = strfun * c_a1 - (c_a2 * epsilon)
    strfun = 1.0 + (1 / (1.0 - c_a1 * epsilon_0a)) * strfun

    # ======================================================================

    # Critical field at zero temperature and current, corrected for strain
    b_c20_eps = b_c20max * strfun

    # Critical temperature at zero field and current, corrected for strain
    temp_c0_eps = temp_c0max * strfun ** (1 / 3)

    # If input temperature is over the strain adjusted critical temperature then report error
    if temp_conductor / temp_c0_eps >= 1.0:
        eh.fdiags[0] = temp_conductor
        eh.fdiags[1] = temp_c0_eps
        eh.report_error(159)

    # Reduced temperature at zero field, corrected for strain
    # f_temp_conductor_critical > 1 is permitted, indicating the temperature is above the critical value at zero field.
    f_temp_conductor_critical_no_field = temp_conductor / temp_c0_eps

    # If input field is over the strain adjusted critical field then report error
    if b_conductor / b_c20_eps >= 1.0:
        eh.fdiags[0] = b_conductor
        eh.fdiags[1] = b_c20_eps
        eh.report_error(160)

    # Reduced field at zero temperature, taking account of strain
    f_b_conductor_critical_no_temp = b_conductor / b_c20_eps

    # Critical temperature at given strain and field
    # tcrit is not used in the calculation of jcrit.
    if (
        f_b_conductor_critical_no_temp < 1.0
    ):  # Normal case, field is within critical surface
        temp_critical = temp_c0_eps * (1.0 - f_b_conductor_critical_no_temp) ** (
            1 / 1.52
        )
    else:  # Abnormal case, field is too high.
        temp_critical = -temp_c0_eps * abs(1.0 - f_b_conductor_critical_no_temp) ** (
            1 / 1.52
        )  # Prevents NaNs. Sets tcrit negative

    # Critical field (T) at given strain and temperature
    b_critical = b_c20_eps * (1.0 - f_temp_conductor_critical_no_field**1.52)

    jc1 = (csc / b_conductor) * strfun

    # Check if we are inside the critical surface
    if (
        (f_temp_conductor_critical_no_field > 0)
        and (f_temp_conductor_critical_no_field < 1)
        and (b_conductor > 0)
        and (b_conductor < b_critical)
        and (b_critical > 0)
    ):
        # Reduced field at given strain and temperature
        b_reduced = b_conductor / b_critical

        jc2 = (1.0 - f_temp_conductor_critical_no_field**1.52) * (
            1.0 - f_temp_conductor_critical_no_field**2
        )
        jc3 = b_reduced**p * (1.0 - b_reduced) ** q
        j_scaling = jc1 * jc2 * jc3

    else:
        # Outside the critical surface.
        # We construct a simple function that is always negative and
        # becomes more negative as field and temperature increase.
        jc2 = f_temp_conductor_critical_no_field
        jc3 = b_conductor / max(b_critical, 1.0e-8)
        j_scaling = -abs(jc1 * jc2 * jc3)

    return j_scaling, b_critical, temp_critical


def croco(j_crit_sc, conductor_area, croco_od, croco_thick):
    """'CroCo' (cross-conductor) strand and cable design for
    'REBCO' 2nd generation HTS superconductor
    Updated 13/11/18 using data from Lewandowska et al 2018.
    """
    d = croco_od
    # d = conductor_width / 3.0d0 - thwcndut * ( 2.0d0 / 3.0d0 )

    croco_id = d - 2.0 * croco_thick  # scaling * 5.4d-3
    if croco_id <= 0.0:
        logger.warning("Negitive inner croco diameter")

    # Define the scaling factor for the input REBCO variable
    # Ratio of new croco inner diameter and fixed base line value
    scaling = croco_id / 5.4e-3
    tape_width = scaling * 3.75e-3
    # Properties of a single strand
    tape_thickness = (
        rebco_variables.rebco_thickness
        + rebco_variables.copper_thick
        + rebco_variables.hastelloy_thickness
    )
    stack_thickness = np.sqrt(croco_id**2 - tape_width**2)
    tapes = stack_thickness / tape_thickness

    copper_area = (
        np.pi * croco_thick * d
        - np.pi * croco_thick**2
        + rebco_variables.copper_thick
        * tape_width
        * tapes  # copper tube  # copper in tape
    )
    hastelloy_area = rebco_variables.hastelloy_thickness * tape_width * tapes
    solder_area = np.pi / 4.0 * croco_id**2 - stack_thickness * tape_width

    rebco_area = rebco_variables.rebco_thickness * tape_width * tapes
    croco_strand_area = np.pi / 4.0 * d**2
    croco_strand_critical_current = j_crit_sc * rebco_area

    # Conductor properties
    # conductor%number_croco = conductor%acs*(1.0-cable_helium_fraction-copper_bar)/croco_strand_area
    conductor_critical_current = croco_strand_critical_current * 6.0
    # Area of core = area of strand
    conductor_copper_bar_area = croco_strand_area
    conductor_copper_area = copper_area * 6.0 + conductor_copper_bar_area
    conductor_copper_fraction = conductor_copper_area / conductor_area

    # Helium area is set by the user.
    # conductor_helium_area = cable_helium_fraction * conductor_acs
    conductor_helium_area = np.pi / 2.0 * d**2
    conductor_helium_fraction = conductor_helium_area / conductor_area

    conductor_hastelloy_area = hastelloy_area * 6.0
    conductor_hastelloy_fraction = conductor_hastelloy_area / conductor_area

    conductor_solder_area = solder_area * 6.0
    conductor_solder_fraction = conductor_solder_area / conductor_area

    conductor_rebco_area = rebco_area * 6.0
    conductor_rebco_fraction = conductor_rebco_area / conductor_area

    return (
        croco_strand_area,
        croco_strand_critical_current,
        conductor_copper_area,
        conductor_copper_fraction,
        conductor_copper_bar_area,
        conductor_hastelloy_area,
        conductor_hastelloy_fraction,
        conductor_helium_area,
        conductor_helium_fraction,
        conductor_solder_area,
        conductor_solder_fraction,
        conductor_rebco_area,
        conductor_rebco_fraction,
        conductor_critical_current,
    )


def current_density_margin(ttest, isumat, jsc, bmax, strain, bc20m, tc0m, c0=None):
    """Current density margin is the difference between the operating current density and
    the critical current density of a superconductor, at a given temperature and field.
    It is zero at the current-sharing temperature.
    ttest : input real :    Temperature
    isumat : input real :   Switch for superconductor material
                            (This has different global names depending on which coil is referred to.)
    jsc : input real : actual current density
    bmax : input real : magnetic field (T)
    strain : input real : superconductor strain
    bc20m, tc0m : input real : superconductor parameters
    """

    # Critical current density jcrit
    if isumat == 1:
        jcrit, _, _ = itersc(ttest, bmax, strain, bc20m, tc0m)
    elif isumat == 3:
        jcrit, _ = jcrit_nbti(ttest, bmax, c0, bc20m, tc0m)
    if isumat == 4:
        jcrit, _, _ = itersc(ttest, bmax, strain, bc20m, tc0m)
    elif isumat == 5:
        jcrit, _, _ = western_superconducting_nb3sn(ttest, bmax, strain, bc20m, tc0m)
    elif isumat == 7:
        jcrit, _, _ = gl_nbti(ttest, bmax, strain, bc20m, tc0m)
    elif isumat == 8:
        jcrit, _, _ = gl_rebco(ttest, bmax, strain, bc20m, tc0m)
    elif isumat == 9:
        jcrit, _, _ = hijc_rebco(ttest, bmax, strain, bc20m, tc0m)

    return jcrit - jsc
