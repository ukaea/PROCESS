import numpy as np

from process.fortran import error_handling as eh


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
