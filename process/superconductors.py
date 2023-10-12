import numpy as np

from process.fortran import error_handling as eh, rebco_variables


def itersc(thelium, bmax, strain, bc20max, tc0max):
    """Implementation of ITER Nb3Sn critical surface implementation
    author: R Kemp, CCFE, Culham Science Centre
    author: P J Knight, CCFE, Culham Science Centre
    thelium : input real : Coolant/SC temperature (K)
    bmax : input real : Magnetic field at conductor (T)
    strain : input real : Strain in superconductor
    bc20max : input real : Upper critical field (T) for superconductor
    at zero temperature and strain
    tc0max : input real : Critical temperature (K) at zero field and strain
    jcrit : output real : Critical current density in superconductor (A/m2)
    bcrit : output real : Critical field (T)
    tcrit : output real : Critical temperature (K)
    This routine calculates the critical current density and
    temperature in the superconducting TF coils using the
    ITER Nb3Sn critical surface model.
    $J_C(B,T,\\epsilon)$ Parameterization for ITER Nb3Sn production,
    L. Bottura, CERN-ITER Collaboration Report, Version 2, April 2nd 2008
    (distributed by Arnaud Devred, ITER, 10th April 2008)
    ITER Nb3Sn critical surface parameterization (2MMF7J) (2008),
    https://user.iter.org/?uid=2MMF7J&action=get_document
    ITER DDD 11-7: Magnets - conductors (2NBKXY) (2009),
    https://user.iter.org/?uid=2NBKXY&action=get_document
    """

    csc = 19922.0e6  # scaling constant C
    p = 0.63  # low field exponent p
    q = 2.1  # high field exponent q
    ca1 = 44.48  # strain fitting constant C_{a1}
    ca2 = 0.0  # strain fitting constant C_{a2}
    eps0a = 0.00256  # epsilon_{0,a}
    diter = 0.82  # ITER strand diameter (mm)
    cuiter = 0.5  # ITER strand copper fraction

    #  $\epsilon_{sh}$
    epssh = (ca2 * eps0a) / (np.sqrt(ca1**2 - ca2**2))

    # Strain function $s(\epsilon)$
    # 0.83 < s < 1.0, for -0.005 < strain < 0.005
    strfun = np.sqrt(epssh**2 + eps0a**2) - np.sqrt(
        (strain - epssh) ** 2 + eps0a**2
    )
    strfun = strfun * ca1 - ca2 * strain
    strfun = 1.0 + (1 / (1.0 - ca1 * eps0a)) * strfun

    # $B^*_{C2} (0,\epsilon)$
    bc20eps = bc20max * strfun

    # $T^*_C (0,\epsilon)$
    tc0eps = tc0max * strfun ** (1 / 3)

    # Reduced temperature, restricted to be < 1
    # Should remain < 1 for thelium < 0.94*tc0max (i.e. 15 kelvin for i_tf_sc_mat=1)

    if thelium / tc0eps >= 1.0:
        eh.fdiags[0] = thelium
        eh.fdiags[1] = tc0eps
        eh.report_error(159)

    t = min(thelium / tc0eps, 0.9999)

    # Reduced magnetic field at zero temperature
    # Should remain < 1 for bmax < 0.83*bc20max (i.e. 27 tesla for i_tf_sc_mat=1)

    if bmax / bc20eps >= 1.0:
        eh.fdiags[0] = bmax
        eh.fdiags[1] = bc20eps
        eh.report_error(160)

    bzero = min(bmax / bc20eps, 0.9999)

    # Critical temperature (K)
    tcrit = tc0eps * (1.0 - bzero) ** (1 / 1.52)  # bzero must be < 1 to avoid NaNs

    # Critical field (T)
    bcrit = bc20eps * (1.0 - t**1.52)

    # Reduced magnetic field, restricted to be < 1
    if bmax / bcrit >= 1.0:
        eh.fdiags[0] = bmax
        eh.fdiags[1] = bcrit
        eh.report_error(161)

    bred = min(bmax / bcrit, 0.9999)

    #  Critical current density in superconductor (A/m2)
    #  ITER parameterization is for the current in a single strand,
    #  not per unit area, so scalefac converts to current density

    scalefac = np.pi * (0.5 * diter) ** 2 * (1.0 - cuiter)

    jc1 = (csc / bmax) * strfun
    jc2 = (1.0 - t**1.52) * (1.0 - t**2)  # t must be < 1 to avoid NaNs
    jc3 = bred**p * (1.0 - bred) ** q  # bred must be < 1 to avoid NaNs

    jcrit = jc1 * jc2 * jc3 / scalefac

    return jcrit, bcrit, tcrit


def jcrit_nbti(temperature, bmax, c0, bc20max, tc0max):
    """Critical current density in a NbTi superconductor strand
    author: P J Knight, CCFE, Culham Science Centre
    temperature : input real : SC temperature (K)
    bmax : input real : Magnetic field at conductor (T)
    c0   : input real : Scaling constant (A/m2)
    bc20max : input real : Upper critical field (T) for superconductor
    at zero temperature and strain
    tc0max : input real : Critical temperature (K) at zero field and strain
    jcrit : output real : Critical current density in superconductor (A/m2)
    tcrit : output real : Critical temperature (K)
    This routine calculates the critical current density and
    temperature in superconducting TF coils using NbTi
    as the superconductor.
    """

    bratio = bmax / bc20max

    if bmax < bc20max:
        #  Critical temperature (K)
        tcrit = tc0max * (1.0 - bratio) ** 0.59
    else:
        # Allow bmax > bc20max but set error flag
        # Fudge to give real (negative) value if bratio < 1
        tcrit = tc0max * (1.0 - bratio)

    # Allow tbar to be negative but set error flag
    tbar = 1.0 - temperature / tcrit

    #  Critical current density (A/m2)
    jcrit = c0 * (1.0 - bratio) * tbar

    # if ((temperature > tcrit).or.(bmax > bc20max))then
    #     write(*,*)'jcrit_nbti: out of range: ', 'bmax =', bmax, ' bc20max =', bc20max, &
    #               ' temperature =',temperature, ' tcrit =',tcrit
    # end if

    return jcrit, tcrit


def bi2212(bmax, jstrand, tsc, fhts):
    """Fitted parameterization to Bi-2212 superconductor properties
    author: P J Knight, CCFE, Culham Science Centre
    author: M Kovari, CCFE, Culham Science Centre
    bmax    : input real : Magnetic field at conductor (T)
    jstrand : input real : Current density in strand (A/m2)
    tsc     : input real : Superconductor temperature (K)
    fhts    : input real : Adjustment factor (<= 1) to account for strain,
    radiation damage, fatigue or AC losses
    jcrit : output real : Critical current density in strand (A/m2)
    tmarg : output real : Temperature margin (K)
    This routine calculates the critical current density and
    the temperature margin for Bi-2212 superconductor in the TF coils
    using a fit by M. Kovari to measurements described in the reference,
    specifically from the points shown in Figure 6.
    <P>Bi-2212 (Bi<SUB>2</SUB>Sr<SUB>2</SUB>CaCu<SUB>2</SUB>O<SUB>8-x</SUB>)
    is a first-generation high temperature superconductor; it still needs
    to be operated below about 10K, but remains superconducting at much
    higher fields at that temperature than Nb3Sn etc.
    The model's range of validity is T &lt; 20K, adjusted field
    b &lt; 104 T, B &gt; 6 T.
    A transformative superconducting magnet technology for fields well
    above 30 T using isotropic round wire multifilament
    Bi2Sr2CaCu2O8-x conductor, D. C. Larbalestier et al., preprint,
    9th April 2013
    """
    b = bmax / np.exp(-0.168 * (tsc - 4.2))

    #  Engineering (i.e. strand) critical current density (A/m2)

    jcrit = fhts * (1.175e9 * np.exp(-0.02115 * b) - 1.288e8)

    #  Temperature margin (K)
    #  Simple inversion of above calculation, using actual current density
    #  in strand instead of jcrit

    tmarg = (
        1.0
        / 0.168
        * np.log(np.log(1.175e9 / (jstrand / fhts + 1.288e8)) / (0.02115 * bmax))
        + 4.2
        - tsc
    )

    #  Check if ranges of validity have been violated

    if (tsc > 20.0) or (bmax < 6.0) or (b > 104.0):
        eh.fdiags[0] = tsc
        eh.fdiags[1] = bmax
        eh.fdiags[2] = b
        eh.report_error(106)

    return jcrit, tmarg


def gl_nbti(thelium, bmax, strain, bc20max, t_c0):
    """Author: S B L Chislett-McDonald Durham University
    Category: subroutine

    Critical current density of the superconductor in an ITER
    Nb-Ti strand based on the Ginzburg-Landau theory of superconductivity

    \\begin{equation}
    J_{c,TS}(B,T,\\epsilon_{I}) = A(\\epsilon_{I}) \\left[T_{c}(\\epsilon_{I})*(1-t^2)\\right]^2\\left
    [B_{c2}(\\epsilon_I)*(1-t^\\nu)\\right]^{n-3}b^{p-1}(1-b)^q~.
    \\end{equation}

    - \\( \\thelium \\) -- Coolant/SC temperature [K]
    - \\( \\bmax \\) -- Magnetic field at conductor [T]
    - \\( \\epsilon_{I} \\) -- Intrinsic strain in superconductor [\\%]
    - \\( \\B_{c2}(\\epsilon_I) \\) -- Strain dependent upper critical field [T]
    - \\( \\b \\) -- Reduced field = bmax / \\B_{c2}(\\epsilon_I)*(1-t^\\nu) [unitless]
    - \\( \\T_{c}(\\epsilon_{I}) \\) -- Strain dependent critical temperature (K)
    - \\( \\t \\) -- Reduced temperature = thelium / \\T_{c}(\\epsilon_{I}) [unitless]
    - \\( \\A(\\epsilon_{I}) \\) -- Strain dependent Prefactor [A / ( m\\(^2\\) K\\(^-2) T\\(^n-3))]
    - \\( \\J_{c,TS} \\) --  Critical current density in superconductor [A / m\\(^-2\\)]
    """

    A_0 = 1102e6
    p = 0.49
    q = 0.56
    n = 1.83
    v = 1.42
    c2 = -0.0025
    c3 = -0.0003
    c4 = -0.0001
    em = -0.002e-2
    u = 0.0
    w = 2.2

    epsilon_I = strain - em

    strain_func = (
        1 + c2 * (epsilon_I) ** 2 + c3 * (epsilon_I) ** 3 + c4 * (epsilon_I) ** 4
    )

    T_e = t_c0 * strain_func ** (1 / w)

    t_reduced = thelium / T_e

    A_e = A_0 * strain_func ** (u / w)

    # Critical Field
    bcrit = bc20max * (1 - t_reduced**v) * strain_func

    b_reduced = bmax / bcrit

    # Critical temperature (K)
    tcrit = T_e

    # Critical current density (A/m2)
    if b_reduced <= 1.0:
        jcrit = (
            A_e
            * (T_e * (1 - t_reduced**2)) ** 2
            * bcrit ** (n - 3)
            * b_reduced ** (p - 1)
            * (1 - b_reduced) ** q
        )
    else:  # Fudge to yield negative single valued function of Jc for fields above Bc2
        jcrit = (
            A_e
            * (T_e * (1 - t_reduced**2)) ** 2
            * bcrit ** (n - 3)
            * b_reduced ** (p - 1)
            * (1 - b_reduced) ** 1.0
        )

    return jcrit, bcrit, tcrit


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
    A_0 = 2.95e2
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

    epsilon_I = strain - em

    strain_func = (
        1 + c2 * (epsilon_I) ** 2 + c3 * (epsilon_I) ** 3 + c4 * (epsilon_I) ** 4
    )

    T_e = t_c0 * strain_func ** (1 / w)

    t_reduced = thelium / T_e

    A_e = A_0 * strain_func ** (u / w)

    #  Critical Field
    bcrit = bc20max * (1 - t_reduced) ** s * strain_func

    b_reduced = bmax / bcrit

    #  Critical temperature (K)
    tcrit = T_e

    #  Critical current density (A/m2)
    jcrit = (
        A_e
        * (T_e * (1 - t_reduced**2)) ** 2
        * bcrit ** (n - 3)
        * b_reduced ** (p - 1)
        * (1 - b_reduced) ** q
    )

    return jcrit, bcrit, tcrit


def hijc_rebco(thelium, bmax, strain, bc20max, t_c0):
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
    June 2018, Art no. 4802404, doi: 10.1109/TASC.2018.2815767. And on the experimental
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
    """

    a = 1.4
    b = 2.005
    # critical current density prefactor
    A_0 = 2.2e8
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
    A_t = A_0 + (u * thelium**2) + (v * thelium)

    # Critical current density (A/m2)

    jcrit = (A_t / bmax) * bcrit**b * (bmax / bcrit) ** p * (1 - bmax / bcrit) ** q

    # Jc times HTS area: default area is width 4mm times HTS layer thickness 1 um,
    # divided by the tape area to provide engineering Jc per tape, then multiplied by fraction 0.4
    # to reach the level of current density expected in the space where the tapes are wound in A/m^2!
    jcrit = (
        jcrit
        * (rebco_variables.tape_width * rebco_variables.rebco_thickness)
        / (rebco_variables.tape_width * rebco_variables.tape_thickness)
        * 0.4
    )

    return jcrit, bcrit, tcrit
