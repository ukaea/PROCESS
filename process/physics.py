import numpy as np
import numba as nb
import scipy.integrate as integrate
from scipy.optimize import root_scalar
import math
import process.physics_functions as physics_funcs
from process.utilities.f2py_string_patch import f2py_compatible_to_string
from process.fortran import (
    constraint_variables,
    reinke_variables,
    reinke_module,
    impurity_radiation_module,
    constants,
    physics_functions_module,
    physics_variables,
    physics_module,
    pulse_variables,
    times_variables,
    current_drive_variables,
    error_handling,
    fwbs_variables,
    build_variables,
    divertor_variables,
    numerics,
    stellarator_variables,
    process_output as po,
    profiles_module,
)


@nb.jit(nopython=True, cache=True)
def coulg(j, tempe, ne):
    """Coulomb logarithm
    author: P J Knight, CCFE, Culham Science Centre
    j  : input integer : radial element index in range 1 to nr
    This function calculates the Coulomb logarithm, valid
    for e-e collisions (T_e > 0.01 keV), and for
    e-i collisions (T_e > 0.01*Zeff^2) (Alexander, 9/5/1994).
    <P>The code was supplied by Emiliano Fable, IPP Garching
    (private communication).
    C. A. Ordonez and M. I. Molina, Phys. Plasmas <B>1</B> (1994) 2515
    Rev. Mod. Phys., V.48, Part 1 (1976) 275
    """
    return 15.9 - 0.5 * np.log(ne[j - 1]) + np.log(tempe[j - 1])


@nb.jit(nopython=True, cache=True)
def nuee(j, tempe, ne):
    """Frequency of electron-electron collisions
    author: P J Knight, CCFE, Culham Science Centre
    j  : input integer : radial element index in range 1 to nr
    This function calculates the frequency of electron-electron
    collisions (Hz): <I>NUEE = 4*SQRT(pi)/3*Ne*e**4*lambd/
    SQRT(Me)/Te**1.5</I>
    <P>The code was supplied by Emiliano Fable, IPP Garching
    (private communication).
    Yushmanov, 25th April 1987 (?),
    updated by Pereverzev, 9th November 1994 (?)
    """
    return (
        670.0 * coulg(j, tempe, ne) * ne[j - 1] / (tempe[j - 1] * np.sqrt(tempe[j - 1]))
    )


@nb.jit(nopython=True, cache=True)
def nui(j, zmain, ni, tempi, amain):
    """Full frequency of ion collisions
    author: P J Knight, CCFE, Culham Science Centre
    j  : input integer : radial element index in range 1 to nr
    This function calculates the full frequency of ion
    collisions (Hz).
    <P>The code was supplied by Emiliano Fable, IPP Garching
    (private communication).
    None
    """
    # Coulomb logarithm = 15 is used
    return (
        zmain[j - 1] ** 4
        * ni[j - 1]
        * 322.0
        / (tempi[j - 1] * np.sqrt(tempi[j - 1] * amain[j - 1]))
    )


@nb.jit(nopython=True, cache=True)
def _plascar_bpol(aspect, eps, kappa, delta):
    # Original coding, only suitable for TARTs [STAR Code]

    c1 = kappa**2 / (1.0 + delta) + delta
    c2 = kappa**2 / (1.0 - delta) - delta

    d1 = (kappa / (1.0 + delta)) ** 2 + 1.0
    d2 = (kappa / (1.0 - delta)) ** 2 + 1.0

    c1_aspect = (c1 * eps - 1.0) if aspect < c1 else (1.0 - c1 * eps)

    y1 = np.sqrt(c1_aspect / (1.0 + eps)) * (1.0 + delta) / kappa
    y2 = np.sqrt((c2 * eps + 1.0) / (1.0 - eps)) * (1.0 - delta) / kappa

    h2 = (1.0 + (c2 - 1.0) * eps / 2.0) / np.sqrt((1.0 - eps) * (c2 * eps + 1.0))
    f2 = (d2 * (1.0 - delta) * eps) / ((1.0 - eps) * (c2 * eps + 1.0))
    g = eps * kappa / (1.0 - eps * delta)
    ff2 = f2 * (g + 2.0 * h2 * np.arctan(y2))

    h1 = (1.0 + (1.0 - c1) * eps / 2.0) / np.sqrt((1.0 + eps) * c1_aspect)
    f1 = (d1 * (1.0 + delta) * eps) / ((1.0 + eps) * (c1 * eps - 1.0))

    if aspect < c1:
        ff1 = f1 * (g - h1 * np.log((1.0 + y1) / (1.0 - y1)))
    else:
        ff1 = -f1 * (-g + 2.0 * h1 * np.arctan(y1))

    return ff1, ff2, d1, d2


@nb.jit(nopython=True, cache=True)
def bpol(icurr, ip, qbar, aspect, eps, bt, kappa, delta, perim, rmu0):
    """Function to calculate poloidal field
    author: J Galambos, FEDC/ORNL
    author: P J Knight, CCFE, Culham Science Centre
    icurr  : input integer : current scaling model to use
    ip     : input real :  plasma current (A)
    qbar   : input real :  edge q-bar
    aspect : input real :  plasma aspect ratio
    bt     : input real :  toroidal field on axis (T)
    kappa  : input real :  plasma elongation
    delta  : input real :  plasma triangularity
    perim  : input real :  plasma perimeter (m)
    This function calculates the poloidal field in Tesla,
    using a simple calculation using Stoke's Law for conventional
    tokamaks, or for TARTs, a scaling from Peng, Galambos and
    Shipe (1992).
    J D Galambos, STAR Code : Spherical Tokamak Analysis and Reactor Code,
    unpublished internal Oak Ridge document
    Y.-K. M. Peng, J. Galambos and P.C. Shipe, 1992,
    Fusion Technology, 21, 1729
    """
    if icurr != 2:
        return rmu0 * ip / perim

    ff1, ff2, _, _ = _plascar_bpol(aspect, eps, kappa, delta)

    return bt * (ff1 + ff2) / (2.0 * np.pi * qbar)


@nb.jit(nopython=True, cache=True)
def plasc(qbar, aspect, eps, rminor, bt, kappa, delta):
    """Function to calculate plasma current (Peng scaling)
    author: J Galambos, FEDC/ORNL
    author: P J Knight, CCFE, Culham Science Centre
    aspect : input real :  plasma aspect ratio
    bt     : input real :  toroidal field on axis (T)
    delta  : input real :  plasma triangularity
    kappa  : input real :  plasma elongation
    qbar   : input real :  edge q-bar
    rminor : input real :  plasma minor radius (m)
    This function calculates the plasma current in MA,
    using a scaling from Peng, Galambos and Shipe (1992).
    It is primarily used for Tight Aspect Ratio Tokamaks and is
    selected via <CODE>icurr=2</CODE>.
    J D Galambos, STAR Code : Spherical Tokamak Analysis and Reactor Code,
    unpublished internal Oak Ridge document
    Y.-K. M. Peng, J. Galambos and P.C. Shipe, 1992,
    Fusion Technology, 21, 1729
    """

    ff1, ff2, d1, d2 = _plascar_bpol(aspect, eps, kappa, delta)

    e1 = 2.0 * kappa / (d1 * (1.0 + delta))
    e2 = 2.0 * kappa / (d2 * (1.0 - delta))

    return (
        rminor
        * bt
        / qbar
        * 5.0
        * kappa
        / (2.0 * np.pi**2)
        * (np.arcsin(e1) / e1 + np.arcsin(e2) / e2)
        * (ff1 + ff2)
    )


@nb.jit(nopython=True, cache=True)
def rether(alphan, alphat, dene, dlamie, te, ti, zeffai):
    """Routine to find the equilibration power between the
    ions and electrons
    author: P J Knight, CCFE, Culham Science Centre
    alphan : input real :  density profile index
    alphat : input real :  temperature profile index
    dene   : input real :  electron density (/m3)
    dlamie : input real :  ion-electron coulomb logarithm
    te     : input real :  electron temperature (keV)
    ti     : input real :  ion temperature (keV)
    zeffai : input real :  mass weighted plasma effective charge
    piepv  : output real : ion/electron equilibration power (MW/m3)
    This routine calculates the equilibration power between the
    ions and electrons.
    Unknown origin
    """
    profie = (1.0 + alphan) ** 2 / (
        (2.0 * alphan - 0.5 * alphat + 1.0) * np.sqrt(1.0 + alphat)
    )
    conie = 2.42165e-41 * dlamie * dene**2 * zeffai * profie

    return conie * (ti - te) / (te**1.5)


@nb.jit(nopython=True, cache=True)
def vscalc(
    csawth,
    eps,
    facoh,
    gamma,
    kappa,
    rmajor,
    rplas,
    plascur,
    t_fusion_ramp,
    tburn,
    rli,
    rmu0,
):
    """Volt-second requirements
    author: P J Knight, CCFE, Culham Science Centre
    csawth : input real :  coefficient for sawteeth effects
    eps    : input real :  inverse aspect ratio
    facoh  : input real :  fraction of plasma current produced inductively
    gamma  : input real :  Ejima coeff for resistive start-up V-s component
    kappa  : input real :  plasma elongation
    plascur: input real :  plasma current (A)
    rli    : input real :  plasma normalised inductivity
    rmajor : input real :  plasma major radius (m)
    rplas  : input real :  plasma resistance (ohm)
    t_fusion_ramp  : input real :  heating time (s)
    tburn  : input real :  burn time (s)
    phiint : output real : internal plasma volt-seconds (Wb)
    rlp    : output real : plasma inductance (H)
    vsbrn  : output real : volt-seconds needed during flat-top (heat+burn) (Wb)
    vsind  : output real : internal and external plasma inductance V-s (Wb)
    vsres  : output real : resistive losses in start-up volt-seconds (Wb)
    vsstt  : output real : total volt-seconds needed (Wb)
    This subroutine calculates the volt-second requirements and some
    other related items.
    """
    # Internal inductance

    rlpint = rmu0 * rmajor * rli / 2.0
    phiint = rlpint * plascur

    # Start-up resistive component
    # Uses ITER formula without the 10 V-s add-on

    vsres = gamma * rmu0 * plascur * rmajor

    # Hirshman, Neilson: Physics of Fluids, 29 (1986) p790
    # fit for external inductance

    aeps = (1.0 + 1.81 * np.sqrt(eps) + 2.05 * eps) * np.log(8.0 / eps) - (
        2.0 + 9.25 * np.sqrt(eps) - 1.21 * eps
    )
    beps = (
        0.73 * np.sqrt(eps) * (1.0 + 2.0 * eps**4 - 6.0 * eps**5 + 3.7 * eps**6)
    )
    rlpext = rmajor * rmu0 * aeps * (1.0 - eps) / (1.0 - eps + beps * kappa)

    rlp = rlpext + rlpint

    # Inductive V-s component

    vsind = rlp * plascur
    vsstt = vsres + vsind

    # Loop voltage during flat-top
    # Include enhancement factor in flattop V-s requirement
    # to account for MHD sawtooth effects.

    vburn = plascur * rplas * facoh * csawth

    # N.B. tburn on first iteration will not be correct
    # if the pulsed reactor option is used, but the value
    # will be correct on subsequent calls.

    vsbrn = vburn * (t_fusion_ramp + tburn)
    vsstt = vsstt + vsbrn

    return phiint, rlp, vsbrn, vsind, vsres, vsstt


@nb.jit(nopython=True, cache=True)
def culblm(bt, dnbeta, plascur, rminor):
    """Beta scaling limit
    author: P J Knight, CCFE, Culham Science Centre
    bt      : input real :  toroidal B-field on plasma axis (T)
    dnbeta  : input real :  Troyon-like g coefficient
    plascur : input real :  plasma current (A)
    rminor  : input real :  plasma minor axis (m)
    betalim : output real : beta limit as defined below
    This subroutine calculates the beta limit, using
    the algorithm documented in AEA FUS 172.
    The limit applies to beta defined with respect to the total B-field.
    Switch iculbl determines which components of beta to include.

    If iculbl = 0, then the limit is applied to the total beta
    If iculbl = 1, then the limit is applied to the thermal beta only
    If iculbl = 2, then the limit is applied to the thermal + neutral beam beta components
    If iculbl = 3, then the limit is applied to the toroidal beta

    The default value for the g coefficient is dnbeta = 3.5
    AEA FUS 172: Physics Assessment for the European Reactor Study
    """

    return 0.01 * dnbeta * (plascur / 1.0e6) / (rminor * bt)


@nb.jit(nopython=True, cache=True)
def conhas(alphaj, alphap, bt, delta95, eps, kappa95, p0, rmu0):
    """Routine to calculate the F coefficient used for scaling the
    plasma current
    author: P J Knight, CCFE, Culham Science Centre
    alphaj   : input real :  current profile index
    alphap   : input real :  pressure profile index
    bt       : input real :  toroidal field on axis (T)
    delta95  : input real :  plasma triangularity 95%
    eps      : input real :  inverse aspect ratio
    kappa95  : input real :  plasma elongation 95%
    p0       : input real :  central plasma pressure (Pa)
    fq       : output real : scaling for edge q from circular
    cross-section cylindrical case
    This routine calculates the F coefficient used for scaling the
    plasma current, using the Connor-Hastie scaling given in
    AEA FUS 172: Physics Assessment for the European Reactor Study
    """

    # Exponent in Connor-Hastie current profile - matching total
    # current gives the following trivial relation
    lamda = alphaj

    # Exponent in Connor-Hastie pressure profile
    nu = alphap

    # Central plasma beta
    beta0 = 2.0 * rmu0 * p0 / (bt**2)

    # Plasma internal inductance
    lamp1 = 1.0 + lamda
    li = lamp1 / lamda * (lamp1 / lamda * np.log(lamp1) - 1.0)

    # T/r in AEA FUS 172
    kap1 = kappa95 + 1.0
    tr = kappa95 * delta95 / kap1**2

    # E/r in AEA FUS 172
    er = (kappa95 - 1.0) / kap1

    # T primed in AEA FUS 172
    tprime = 2.0 * tr * lamp1 / (1.0 + 0.5 * lamda)

    # E primed in AEA FUS 172
    eprime = er * lamp1 / (1.0 + lamda / 3.0)

    # Delta primed in AEA FUS 172
    deltap = (0.5 * kap1 * eps * 0.5 * li) + (
        beta0 / (0.5 * kap1 * eps)
    ) * lamp1**2 / (1.0 + nu)

    # Delta/R0 in AEA FUS 172
    deltar = beta0 / 6.0 * (1.0 + 5.0 * lamda / 6.0 + 0.25 * lamda**2) + (
        0.5 * kap1 * eps
    ) ** 2 * 0.125 * (1.0 - (lamda**2) / 3.0)

    # F coefficient
    return (0.5 * kap1) ** 2 * (
        1.0
        + eps**2 * (0.5 * kap1) ** 2
        + 0.5 * deltap**2
        + 2.0 * deltar
        + 0.5 * (eprime**2 + er**2)
        + 0.5 * (tprime**2 + 4.0 * tr**2)
    )


def bsinteg(
    y,
    dene,
    ten,
    bt,
    rminor,
    rmajor,
    zeff,
    alphat,
    alphan,
    q0,
    q95,
    betat,
    echarge,
    rmu0,
):
    """Integrand function for Nevins et al bootstrap current scaling
    author: P J Knight, CCFE, Culham Science Centre
    y : input real : abscissa of integration, = normalised minor radius
    This function calculates the integrand function for the
    Nevins et al bootstrap current scaling, 4/11/90.
    """
    # Constants for fit to q-profile
    c1 = 1.0
    c2 = 1.0
    c3 = 1.0

    # Compute average electron beta
    betae = dene * ten * 1.0e3 * echarge / (bt**2 / (2.0 * rmu0))

    nabla = rminor * np.sqrt(y) / rmajor
    x = (1.46 * np.sqrt(nabla) + 2.4 * nabla) / (1.0 - nabla) ** 1.5
    z = zeff
    d = (
        1.414 * z
        + z * z
        + x * (0.754 + 2.657 * z + 2.0 * z * z)
        + x * x * (0.348 + 1.243 * z + z * z)
    )
    al2 = -x * (0.884 + 2.074 * z) / d
    a2 = alphat * (1.0 - y) ** (alphan + alphat - 1.0)
    alphai = -1.172 / (1.0 + 0.462 * x)
    a1 = (alphan + alphat) * (1.0 - y) ** (alphan + alphat - 1.0)
    al1 = x * (0.754 + 2.21 * z + z * z + x * (0.348 + 1.243 * z + z * z)) / d

    # q-profile
    q = q0 + (q95 - q0) * (c1 * y + c2 * y * y + c3 * y**3) / (c1 + c2 + c3)

    pratio = (betat - betae) / betae

    return (q / q95) * (al1 * (a1 + pratio * (a1 + alphai * a2)) + al2 * a2)


def diamagnetic_fraction_hender(beta):
    """author: S.I. Muldrew, CCFE, Culham Science Centre
    Diamagnetic contribution at tight aspect ratio.
    Tim Hender fit
    """
    return beta / 2.8


def diamagnetic_fraction_scene(beta, q95, q0):
    """author: S.I. Muldrew, CCFE, Culham Science Centre
    Diamagnetic fraction based on SCENE fit by Tim Hender
    See Issue #992
    """
    return beta * (0.1 * q95 / q0 + 0.44) * 4.14e-1


def ps_fraction_scene(beta):
    """author: S.I. Muldrew, CCFE, Culham Science Centre
    Pfirsch-Schlüter fraction based on SCENE fit by Tim Hender
    See Issue #992
    """
    return -9e-2 * beta


@nb.jit(nopython=True, cache=True)
def nuis(j, rmajor, mu, sqeps, tempi, amain, zmain, ni):
    """Relative frequency of ion collisions
    author: P J Knight, CCFE, Culham Science Centre
    j  : input integer : radial element index in range 1 to nr
    This function calculates the relative frequency of ion
    collisions: <I>NU* = Nui*q*Rt/eps**1.5/Vti</I>
    The full ion collision frequency NUI is used.
    <P>The code was supplied by Emiliano Fable, IPP Garching
    (private communication).
    Yushmanov, 30th April 1987 (?)
    """
    return (
        3.2e-6
        * nui(j, zmain, ni, tempi, amain)
        * rmajor
        / (
            np.abs(mu[j - 1] + 1.0e-4)
            * sqeps[j - 1] ** 3
            * np.sqrt(tempi[j - 1] / amain[j - 1])
        )
    )


@nb.jit(nopython=True, cache=True)
def dcsa(j, nr, rmajor, bt, triang, ne, ni, tempe, tempi, mu, rho, zef, sqeps):
    """Grad(ln(ne)) coefficient in the Sauter bootstrap scaling
    author: P J Knight, CCFE, Culham Science Centre
    j  : input integer : radial element index in range 2 to nr
    nr : input integer : maximum value of j
    This function calculates the coefficient scaling grad(ln(ne))
    in the Sauter bootstrap current scaling.
    Code by Angioni, 29th May 2002.
    <P>The code was supplied by Emiliano Fable, IPP Garching
    (private communication).
    O. Sauter, C. Angioni and Y. R. Lin-Liu,
    Physics of Plasmas <B>6</B> (1999) 2834
    O. Sauter, C. Angioni and Y. R. Lin-Liu, (ERRATA)
    Physics of Plasmas <B>9</B> (2002) 5140

    DCSA $\\equiv \\mathcal{L}_{31}$, Eq.14a, Sauter et al, 1999
    """
    zz = zef[j - 1]
    zft = tpf(j, triang, sqeps)
    _nues = nues(j, rmajor, zef, mu, sqeps, tempe, ne)
    zdf = (1.0 + (1.0 - 0.1 * zft) * np.sqrt(_nues)) + (0.5 * (1.0 - zft) * _nues) / zz
    zft /= zdf  # $f^{31}_{teff}(\nu_{e*})$, Eq.14b
    zft2 = zft**2
    zft3 = zft2 * zft
    dcsa = (
        ((1.0 + 1.4 / (zz + 1.0)) * zft)
        - (1.9 / (zz + 1.0) * zft2)
        + ((0.3 * zft3 + 0.2 * zft3 * zft) / (zz + 1.0))
    )

    # Corrections suggested by Fable, 15/05/2015
    return dcsa * beta_poloidal_local_total(
        j, nr, rmajor, bt, ne, ni, tempe, tempi, mu, rho
    )


@nb.jit(nopython=True, cache=True)
def hcsa(j, nr, rmajor, bt, triang, ne, ni, tempe, tempi, mu, rho, zef, sqeps):
    """Grad(ln(Te)) coefficient in the Sauter bootstrap scaling
    author: P J Knight, CCFE, Culham Science Centre
    j  : input integer : radial element index in range 2 to nr
    nr : input integer : maximum value of j
    This function calculates the coefficient scaling grad(ln(Te))
    in the Sauter bootstrap current scaling.
    Code by Angioni, 29th May 2002.
    <P>The code was supplied by Emiliano Fable, IPP Garching
    (private communication).
    O. Sauter, C. Angioni and Y. R. Lin-Liu,
    Physics of Plasmas <B>6</B> (1999) 2834
    O. Sauter, C. Angioni and Y. R. Lin-Liu, (ERRATA)
    Physics of Plasmas <B>9</B> (2002) 5140
    """
    zz = zef[j - 1]
    zft = tpf(j, triang, sqeps)
    _nues = nues(j, rmajor, zef, mu, sqeps, tempe, ne)
    _nues_sqrt = np.sqrt(_nues)
    zdf = (1.0 + 0.26 * (1.0 - zft) * _nues_sqrt) + (
        0.18 * (1.0 - 0.37 * zft) * _nues / np.sqrt(zz)
    )
    zfte = zft / zdf  # $f^{32\_ee}_{teff}(\nu_{e*})$, Eq.15d
    zfte2 = zfte * zfte
    zfte3 = zfte * zfte2
    zfte4 = zfte2 * zfte2

    zdf = (1.0 + (1.0 + 0.6 * zft) * _nues_sqrt) + (
        0.85 * (1.0 - 0.37 * zft) * _nues * (1.0 + zz)
    )
    zfti = zft / zdf  # $f^{32\_ei}_{teff}(\nu_{e*})$, Eq.15e
    zfti2 = zfti * zfti
    zfti3 = zfti * zfti2
    zfti4 = zfti2 * zfti2

    # $F_{32\_ee}(X)$, Eq.15b
    hcee = (
        ((0.05 + 0.62 * zz) / zz / (1.0 + 0.44 * zz) * (zfte - zfte4))
        + ((zfte2 - zfte4 - 1.2 * (zfte3 - zfte4)) / (1.0 + 0.22 * zz))
        + (1.2 / (1.0 + 0.5 * zz) * zfte4)
    )

    # $F_{32\_ei}(Y)$, Eq.15c
    hcei = (
        (-(0.56 + 1.93 * zz) / zz / (1.0 + 0.44 * zz) * (zfti - zfti4))
        + (4.95 / (1.0 + 2.48 * zz) * (zfti2 - zfti4 - 0.55 * (zfti3 - zfti4)))
        - (1.2 / (1.0 + 0.5 * zz) * zfti4)
    )

    # Corrections suggested by Fable, 15/05/2015
    return beta_poloidal_local(j, nr, rmajor, bt, ne, tempe, mu, rho) * (
        hcee + hcei
    ) + dcsa(
        j, nr, rmajor, bt, triang, ne, ni, tempe, tempi, mu, rho, zef, sqeps
    ) * beta_poloidal_local(
        j, nr, rmajor, bt, ne, tempe, mu, rho
    ) / beta_poloidal_local_total(
        j, nr, rmajor, bt, ne, ni, tempe, tempi, mu, rho
    )


@nb.jit(nopython=True, cache=True)
def xcsa(
    j,
    nr,
    rmajor,
    bt,
    triang,
    mu,
    sqeps,
    tempi,
    tempe,
    amain,
    zmain,
    ni,
    ne,
    rho,
    zef,
):
    """Grad(ln(Ti)) coefficient in the Sauter bootstrap scaling
    author: P J Knight, CCFE, Culham Science Centre
    j  : input integer : radial element index in range 2 to nr
    nr : input integer : maximum value of j
    This function calculates the coefficient scaling grad(ln(Ti))
    in the Sauter bootstrap current scaling.
    Code by Angioni, 29th May 2002.
    <P>The code was supplied by Emiliano Fable, IPP Garching
    (private communication).
    O. Sauter, C. Angioni and Y. R. Lin-Liu,
    Physics of Plasmas <B>6</B> (1999) 2834
    O. Sauter, C. Angioni and Y. R. Lin-Liu, (ERRATA)
    Physics of Plasmas <B>9</B> (2002) 5140
    """
    zz = zef[j - 1]
    zft = tpf(j, triang, sqeps)
    _nues = nues(j, rmajor, zef, mu, sqeps, tempe, ne)
    zdf = (1.0 + (1.0 - 0.1 * zft) * np.sqrt(_nues)) + 0.5 * (
        1.0 - 0.5 * zft
    ) * _nues / zz
    zfte = zft / zdf  # $f^{34}_{teff}(\nu_{e*})$, Eq.16b

    # Eq.16a
    xcsa = ((1.0 + 1.4 / (zz + 1.0)) * zfte - 1.9 / (zz + 1.0) * zfte * zfte) + (
        0.3 * zfte * zfte + 0.2 * zfte * zfte * zfte
    ) * zfte / (zz + 1.0)

    # $\alpha_0$, Eq.17a
    a0 = (-1.17 * (1.0 - zft)) / (1.0 - 0.22 * zft - 0.19 * zft * zft)

    _nuis = nuis(j, rmajor, mu, sqeps, tempi, amain, zmain, ni)
    _nuis_sqrt = np.sqrt(_nuis)
    a1 = _nuis**2 * zft**6

    # $\alpha(\nu_{i*})$, Eq.17b
    alp = (
        (a0 + 0.25 * (1.0 - zft * zft) * _nuis_sqrt) / (1.0 + 0.5 * _nuis_sqrt)
        + 0.315 * a1
    ) / (1.0 + 0.15 * a1)

    # Corrections suggested by Fable, 15/05/2015
    return (
        beta_poloidal_local_total(j, nr, rmajor, bt, ne, ni, tempe, tempi, mu, rho)
        - beta_poloidal_local(j, nr, rmajor, bt, ne, tempe, mu, rho)
    ) * (xcsa * alp) + dcsa(
        j, nr, rmajor, bt, triang, ne, ni, tempe, tempi, mu, rho, zef, sqeps
    ) * (
        1.0
        - beta_poloidal_local(j, nr, rmajor, bt, ne, tempe, mu, rho)
        / beta_poloidal_local_total(j, nr, rmajor, bt, ne, ni, tempe, tempi, mu, rho)
    )


@nb.jit(nopython=True, cache=True)
def nues(j, rmajor, zef, mu, sqeps, tempe, ne):
    """Relative frequency of electron collisions
    author: P J Knight, CCFE, Culham Science Centre
    j  : input integer : radial element index in range 1 to nr
    This function calculates the relative frequency of electron
    collisions: <I>NU* = Nuei*q*Rt/eps**1.5/Vte</I>
    The electron-ion collision frequency NUEI=NUEE*1.4*ZEF is
    used.
    <P>The code was supplied by Emiliano Fable, IPP Garching
    (private communication).
    Yushmanov, 30th April 1987 (?)
    """
    return (
        nuee(j, tempe, ne)
        * 1.4
        * zef[j - 1]
        * rmajor
        / np.abs(mu[j - 1] * (sqeps[j - 1] ** 3) * np.sqrt(tempe[j - 1]) * 1.875e7)
    )


@nb.jit(nopython=True, cache=True)
def beta_poloidal_local(j, nr, rmajor, bt, ne, tempe, mu, rho):
    """Local beta poloidal calculation
    author: P J Knight, CCFE, Culham Science Centre
    j  : input integer : radial element index in range 1 to nr
    nr : input integer : maximum value of j
    This function calculates the local beta poloidal.
    <P>The code was supplied by Emiliano Fable, IPP Garching
    (private communication).
    <P>beta poloidal = 4*pi*ne*Te/Bpo**2
    Pereverzev, 25th April 1989 (?)
    """
    return (
        np.where(
            j != nr,
            1.6e-4 * np.pi * (ne[j] + ne[j - 1]) * (tempe[j] + tempe[j - 1]),
            6.4e-4 * np.pi * ne[j - 1] * tempe[j - 1],
        )
        * (rmajor / (bt * rho[j - 1] * np.abs(mu[j - 1] + 1.0e-4))) ** 2
    )


@nb.jit(nopython=True, cache=True)
def beta_poloidal_local_total(j, nr, rmajor, bt, ne, ni, tempe, tempi, mu, rho):
    """Local beta poloidal calculation, including ion pressure
    author: P J Knight, CCFE, Culham Science Centre
    j  : input integer : radial element index in range 1 to nr
    nr : input integer : maximum value of j
    This function calculates the local total beta poloidal.
    <P>The code was supplied by Emiliano Fable, IPP Garching
    (private communication).
    <P>beta poloidal = 4*pi*(ne*Te+ni*Ti)/Bpo**2
    where ni is the sum of all ion densities (thermal)
    Pereverzev, 25th April 1989 (?)
    E Fable, private communication, 15th May 2014
    """

    return (
        np.where(
            j != nr,
            1.6e-4
            * np.pi
            * (
                ((ne[j] + ne[j - 1]) * (tempe[j] + tempe[j - 1]))
                + ((ni[j] + ni[j - 1]) * (tempi[j] + tempi[j - 1]))
            ),
            6.4e-4 * np.pi * (ne[j - 1] * tempe[j - 1] + ni[j - 1] * tempi[j - 1]),
        )
        * (rmajor / (bt * rho[j - 1] * np.abs(mu[j - 1] + 1.0e-4))) ** 2
    )


@nb.jit(nopython=True, cache=True)
def tpf(j, triang, sqeps, fit=1):
    """Trapped particle fraction
    author: P J Knight, CCFE, Culham Science Centre
    j  : input integer : radial element index in range 1 to nr
    fit : input integer : (1)=ASTRA method, 2=Equation from Sauter2002, 3=Equation from Sauter2013
    This function calculates the trapped particle fraction at
    a given radius.
    <P>A number of different fits are provided, but the one
    to be used is hardwired prior to run-time.
    <P>The code was supplied by Emiliano Fable, IPP Garching
    (private communication).
    O. Sauter et al, Plasma Phys. Contr. Fusion <B>44</B> (2002) 1999
    O. Sauter, 2013:
    http://infoscience.epfl.ch/record/187521/files/lrp_012013.pdf
    """
    s = sqeps[j - 1]
    eps = s * s

    if fit == 1:
        # ASTRA method, from Emiliano Fable, private communication
        # (Excluding h term which dominates for inverse aspect ratios < 0.5,
        # and tends to take the trapped particle fraction to 1)

        zz = 1.0 - eps
        return 1.0 - zz * np.sqrt(zz) / (1.0 + 1.46 * s)
    elif fit == 2:
        # Equation 4 of Sauter 2002
        # Similar to, but not quite identical to g above

        return 1.0 - (1.0 - eps) ** 2 / (1.0 + 1.46 * s) / np.sqrt(1.0 - eps * eps)
    elif fit == 3:
        # Includes correction for triangularity

        epseff = 0.67 * (1.0 - 1.4 * triang * np.abs(triang)) * eps

        return 1.0 - np.sqrt((1.0 - eps) / (1.0 + eps)) * (1.0 - epseff) / (
            1.0 + 2.0 * np.sqrt(epseff)
        )

    raise RuntimeError(f"fit={fit} is not valid. Must be 1, 2, or 3")


class Physics:
    def __init__(self, plasma_profile, current_drive):
        self.outfile = constants.nout
        self.plasma_profile = plasma_profile
        self.current_drive = current_drive

    def physics(self):
        """
        Routine to calculate tokamak plasma physics information
        author: P J Knight, CCFE, Culham Science Centre
        None
        This routine calculates all the primary plasma physics
        M. Kovari et al, 2014, "PROCESS": A systems code for fusion power plants -
        Part 1: Physics https://www.sciencedirect.com/science/article/pii/S0920379614005961
        H. Zohm et al, 2013, On the Physics Guidelines for a Tokamak DEMO
        https://iopscience.iop.org/article/10.1088/0029-5515/53/7/073019
        T. Hartmann, 2013, Development of a modular systems code to analyse the
        implications of physics assumptions on the design of a demonstration fusion power plant
        https://inis.iaea.org/search/search.aspx?orig_q=RN:45031642
        """
        # kappaa_IPB = physics_variables.vol / (
        #    2.0e0
        #    * numpy.pi
        #    * numpy.pi
        #    * physics_variables.rminor
        #    * physics_variables.rminor
        #    * physics_variables.rmajor
        # )

        if physics_variables.icurr == 2:
            physics_variables.q95 = (
                physics_variables.q * 1.3e0 * (1.0e0 - physics_variables.eps) ** 0.6e0
            )
        else:
            physics_variables.q95 = (
                physics_variables.q
            )  # i.e. input (or iteration variable) value

        # Calculate plasma composition
        # Issue #261 Remove old radiation model (imprad_model=0)
        self.plasma_composition()

        # Calculate plasma current
        (
            physics_variables.alphaj,
            physics_variables.rli,
            physics_variables.bp,
            physics_variables.qstar,
            physics_variables.plascur,
        ) = self.culcur(
            physics_variables.alphaj,
            physics_variables.alphap,
            physics_variables.bt,
            physics_variables.eps,
            physics_variables.icurr,
            physics_variables.iprofile,
            physics_variables.kappa,
            physics_variables.kappa95,
            physics_variables.p0,
            physics_variables.pperim,
            physics_variables.q0,
            physics_variables.q,
            physics_variables.rli,
            physics_variables.rmajor,
            physics_variables.rminor,
            physics_variables.sf,
            physics_variables.triang,
            physics_variables.triang95,
        )

        # Calculate density and temperature profile quantities
        # If physics_variables.ipedestal = 1 then set pedestal density to
        #   physics_variables.fgwped * Greenwald density limit
        # Note: this used to be done before plasma current
        if (physics_variables.ipedestal == 1) and (physics_variables.fgwped >= 0e0):
            physics_variables.neped = (
                physics_variables.fgwped
                * 1.0e14
                * physics_variables.plascur
                / (np.pi * physics_variables.rminor * physics_variables.rminor)
            )

        if (physics_variables.ipedestal == 1) and (physics_variables.fgwsep >= 0e0):
            physics_variables.nesep = (
                physics_variables.fgwsep
                * 1.0e14
                * physics_variables.plascur
                / (np.pi * physics_variables.rminor * physics_variables.rminor)
            )

        self.plasma_profile.run()

        # Calculate total magnetic field [T]
        physics_variables.btot = physics_functions_module.total_mag_field()

        # Calculate physics_variables.beta poloidal [-]
        physics_variables.betap = physics_functions_module.beta_poloidal()

        # Set PF coil ramp times
        if pulse_variables.lpulse != 1:
            if times_variables.tohsin == 0.0e0:
                times_variables.tohs = physics_variables.plascur / 5.0e5
                times_variables.tramp = times_variables.tohs
                times_variables.tqnch = times_variables.tohs
            else:
                times_variables.tohs = times_variables.tohsin

        else:
            if times_variables.pulsetimings == 0.0e0:
                # times_variables.tramp is input
                times_variables.tohs = physics_variables.plascur / 1.0e5
                times_variables.tqnch = times_variables.tohs

            else:
                # times_variables.tohs is set either in INITIAL or INPUT, or by being
                # iterated using limit equation 41.
                times_variables.tramp = max(times_variables.tramp, times_variables.tohs)
                # tqnch = max(tqnch,tohs)
                times_variables.tqnch = times_variables.tohs

        # Reset second times_variables.tburn value (times_variables.tburn0).
        # This is used to ensure that the burn time is used consistently;
        # see convergence loop in fcnvmc1, evaluators.f90
        times_variables.tburn0 = times_variables.tburn

        # Pulse and down times : The reactor is assumed to be 'down'
        # at all times outside of the plasma current flat-top period.
        # The pulse length is the duration of non-zero plasma current
        times_variables.tpulse = (
            times_variables.tohs
            + times_variables.t_fusion_ramp
            + times_variables.tburn
            + times_variables.tqnch
        )
        times_variables.tdown = (
            times_variables.tramp
            + times_variables.tohs
            + times_variables.tqnch
            + times_variables.tdwell
        )

        # Total cycle time
        times_variables.tcycle = (
            times_variables.tramp
            + times_variables.tohs
            + times_variables.t_fusion_ramp
            + times_variables.tburn
            + times_variables.tqnch
            + times_variables.tdwell
        )

        # Calculate bootstrap current fraction using various models
        current_drive_variables.bscf_iter89 = self.bootstrap_fraction_iter89(
            physics_variables.aspect,
            physics_variables.beta,
            physics_variables.btot,
            current_drive_variables.cboot,
            physics_variables.plascur,
            physics_variables.q95,
            physics_variables.q0,
            physics_variables.rmajor,
            physics_variables.vol,
        )

        betat = (
            physics_variables.beta
            * physics_variables.btot**2
            / physics_variables.bt**2
        )
        current_drive_variables.bscf_nevins = (
            current_drive_variables.cboot
            * self.bootstrap_fraction_nevins(
                physics_variables.alphan,
                physics_variables.alphat,
                betat,
                physics_variables.bt,
                physics_variables.dene,
                physics_variables.plascur,
                physics_variables.q95,
                physics_variables.q0,
                physics_variables.rmajor,
                physics_variables.rminor,
                physics_variables.ten,
                physics_variables.zeff,
            )
        )

        # Wilson scaling uses thermal poloidal beta, not total
        betpth = (
            physics_variables.beta - physics_variables.betaft - physics_variables.betanb
        ) * (physics_variables.btot / physics_variables.bp) ** 2
        current_drive_variables.bscf_wilson = (
            current_drive_variables.cboot
            * self.bootstrap_fraction_wilson(
                physics_variables.alphaj,
                physics_variables.alphap,
                physics_variables.alphat,
                betpth,
                physics_variables.q0,
                physics_variables.q95,
                physics_variables.rmajor,
                physics_variables.rminor,
            )
        )

        # Hender scaling for diamagnetic current at tight physics_variables.aspect ratio
        current_drive_variables.diacf_hender = diamagnetic_fraction_hender(
            physics_variables.beta
        )

        # SCENE scaling for diamagnetic current
        current_drive_variables.diacf_scene = diamagnetic_fraction_scene(
            physics_variables.beta, physics_variables.q95, physics_variables.q0
        )

        # Pfirsch-Schlüter scaling for diamagnetic current
        current_drive_variables.pscf_scene = ps_fraction_scene(physics_variables.beta)

        current_drive_variables.bscf_sauter = (
            current_drive_variables.cboot * self.bootstrap_fraction_sauter()
        )

        current_drive_variables.bscf_sakai = (
            current_drive_variables.cboot
            * self.bootstrap_fraction_sakai(
                betap=physics_variables.betap,
                q95=physics_variables.q95,
                q0=physics_variables.q0,
                alphan=physics_variables.alphan,
                alphat=physics_variables.alphat,
                eps=physics_variables.eps,
                rli=physics_variables.rli,
            )
        )

        if current_drive_variables.bscfmax < 0.0e0:
            current_drive_variables.bootipf = abs(current_drive_variables.bscfmax)
            current_drive_variables.plasipf = current_drive_variables.bootipf
        else:
            if physics_variables.ibss == 1:
                current_drive_variables.bootipf = current_drive_variables.bscf_iter89
            elif physics_variables.ibss == 2:
                current_drive_variables.bootipf = current_drive_variables.bscf_nevins
            elif physics_variables.ibss == 3:
                current_drive_variables.bootipf = current_drive_variables.bscf_wilson
            elif physics_variables.ibss == 4:
                current_drive_variables.bootipf = current_drive_variables.bscf_sauter
            elif physics_variables.ibss == 5:
                # Sakai states that the ACCOME dataset used has the toridal diamagnetic current included in the bootstrap current
                # So the diamagnetic current calculation should be turned off when using, (idia = 0).
                current_drive_variables.bootipf = current_drive_variables.bscf_sakai
            else:
                error_handling.idiags[0] = physics_variables.ibss
                error_handling.report_error(75)

            physics_module.err242 = 0
            if current_drive_variables.bootipf > current_drive_variables.bscfmax:
                current_drive_variables.bootipf = min(
                    current_drive_variables.bootipf, current_drive_variables.bscfmax
                )
                physics_module.err242 = 1

            if physics_variables.idia == 1:
                current_drive_variables.diaipf = current_drive_variables.diacf_hender
            elif physics_variables.idia == 2:
                current_drive_variables.diaipf = current_drive_variables.diacf_scene

            if physics_variables.ips == 1:
                current_drive_variables.psipf = current_drive_variables.pscf_scene

            current_drive_variables.plasipf = (
                current_drive_variables.bootipf
                + current_drive_variables.diaipf
                + current_drive_variables.psipf
            )

        # Plasma driven current fraction (Bootstrap + Diamagnetic
        # + Pfirsch-Schlüter) constrained to be less than
        # or equal to the total fraction of the plasma current
        # produced by non-inductive means (which also includes
        # the current drive proportion)
        physics_module.err243 = 0
        if current_drive_variables.plasipf > physics_variables.fvsbrnni:
            current_drive_variables.plasipf = min(
                current_drive_variables.plasipf, physics_variables.fvsbrnni
            )
            physics_module.err243 = 1

        # Fraction of plasma current produced by inductive means
        physics_variables.facoh = max(1.0e-10, (1.0e0 - physics_variables.fvsbrnni))
        #  Fraction of plasma current produced by auxiliary current drive
        physics_variables.faccd = (
            physics_variables.fvsbrnni - current_drive_variables.plasipf
        )

        # Auxiliary current drive power calculations

        if current_drive_variables.irfcd != 0:
            self.current_drive.cudriv(False)

        # Calculate fusion power + components

        # physics_funcs.palph(self.plasma_profile)
        fusion_rate = physics_funcs.FusionReactionRate(self.plasma_profile)
        fusion_rate.calculate_fusion_rates()
        fusion_rate.set_physics_variables()

        #
        physics_variables.pdt = physics_module.pdtpv * physics_variables.vol
        physics_variables.pdhe3 = physics_module.pdhe3pv * physics_variables.vol
        physics_variables.pdd = physics_module.pddpv * physics_variables.vol

        # Calculate neutral beam slowing down effects
        # If ignited, then ignore beam fusion effects

        if (current_drive_variables.cnbeam != 0.0e0) and (
            physics_variables.ignite == 0
        ):
            (
                physics_variables.betanb,
                physics_variables.dnbeam2,
                physics_variables.palpnb,
            ) = physics_functions_module.beamfus(
                physics_variables.beamfus0,
                physics_variables.betbm0,
                physics_variables.bp,
                physics_variables.bt,
                current_drive_variables.cnbeam,
                physics_variables.dene,
                physics_variables.deni,
                physics_variables.dlamie,
                physics_variables.ealphadt,
                current_drive_variables.enbeam,
                physics_variables.fdeut,
                physics_variables.ftrit,
                current_drive_variables.ftritbm,
                physics_module.sigvdt,
                physics_variables.ten,
                physics_variables.tin,
                physics_variables.vol,
                physics_variables.zeffai,
            )
            physics_variables.fusionrate = (
                physics_variables.fusionrate
                + 1.0e6
                * physics_variables.palpnb
                / (1.0e3 * physics_variables.ealphadt * constants.echarge)
                / physics_variables.vol
            )
            physics_variables.alpharate = (
                physics_variables.alpharate
                + 1.0e6
                * physics_variables.palpnb
                / (1.0e3 * physics_variables.ealphadt * constants.echarge)
                / physics_variables.vol
            )

        physics_variables.pdt = physics_variables.pdt + 5.0e0 * physics_variables.palpnb

        # Create some derived values and add beam contribution to fusion power
        (
            physics_variables.palpmw,
            physics_variables.pneutmw,
            physics_variables.pchargemw,
            physics_variables.betaft,
            physics_variables.palpipv,
            physics_variables.palpepv,
            physics_variables.pfuscmw,
            physics_variables.powfmw,
        ) = physics_functions_module.palph2(
            physics_variables.bt,
            physics_variables.bp,
            physics_variables.dene,
            physics_variables.deni,
            physics_variables.dnitot,
            physics_variables.falpe,
            physics_variables.falpi,
            physics_variables.palpnb,
            physics_variables.ifalphap,
            physics_variables.pchargepv,
            physics_variables.pneutpv,
            physics_variables.ten,
            physics_variables.tin,
            physics_variables.vol,
            physics_variables.palppv,
        )

        # Nominal mean neutron wall load on entire first wall area including divertor and beam holes
        # Note that 'fwarea' excludes these, so they have been added back in.
        if physics_variables.iwalld == 1:
            physics_variables.wallmw = (
                physics_variables.ffwal
                * physics_variables.pneutmw
                / physics_variables.sarea
            )
        else:
            if physics_variables.idivrt == 2:
                # Double null configuration
                physics_variables.wallmw = (
                    (1.0e0 - fwbs_variables.fhcd - 2.0e0 * fwbs_variables.fdiv)
                    * physics_variables.pneutmw
                    / build_variables.fwarea
                )
            else:
                # Single null Configuration
                physics_variables.wallmw = (
                    (1.0e0 - fwbs_variables.fhcd - fwbs_variables.fdiv)
                    * physics_variables.pneutmw
                    / build_variables.fwarea
                )

        # Calculate ion/electron equilibration power

        physics_variables.piepv = rether(
            physics_variables.alphan,
            physics_variables.alphat,
            physics_variables.dene,
            physics_variables.dlamie,
            physics_variables.te,
            physics_variables.ti,
            physics_variables.zeffai,
        )

        # Calculate radiation power

        radpwrdata = physics_funcs.radpwr(self.plasma_profile)
        physics_variables.psyncpv = radpwrdata.psyncpv
        physics_variables.pcoreradpv = radpwrdata.pcoreradpv
        physics_variables.pedgeradpv = radpwrdata.pedgeradpv
        physics_variables.pradpv = radpwrdata.pradpv

        physics_variables.pinnerzoneradmw = (
            physics_variables.pcoreradpv * physics_variables.vol
        )
        physics_variables.pouterzoneradmw = (
            physics_variables.pedgeradpv * physics_variables.vol
        )
        physics_variables.pradmw = physics_variables.pradpv * physics_variables.vol

        # Calculate ohmic power
        (
            physics_variables.pohmpv,
            physics_variables.pohmmw,
            physics_variables.rpfac,
            physics_variables.rplas,
        ) = self.pohm(
            physics_variables.facoh,
            physics_variables.kappa95,
            physics_variables.plascur,
            physics_variables.rmajor,
            physics_variables.rminor,
            physics_variables.ten,
            physics_variables.vol,
            physics_variables.zeff,
        )

        # Calculate L- to H-mode power threshold for different scalings

        physics_variables.pthrmw = physics_functions_module.pthresh(
            physics_variables.dene,
            physics_variables.dnla,
            physics_variables.bt,
            physics_variables.rmajor,
            physics_variables.kappa,
            physics_variables.sarea,
            physics_variables.aion,
            physics_variables.aspect,
        )

        # Enforced L-H power threshold value (if constraint 15 is turned on)

        physics_variables.plhthresh = physics_variables.pthrmw[
            physics_variables.ilhthresh - 1
        ]

        # Power transported to the divertor by charged particles,
        # i.e. excludes neutrons and radiation, and also NBI orbit loss power,
        # which is assumed to be absorbed by the first wall
        if physics_variables.ignite == 0:
            pinj = current_drive_variables.pinjmw
        else:
            pinj = 0.0e0

        physics_variables.pdivt = (
            physics_variables.falpha * physics_variables.palpmw
            + physics_variables.pchargemw
            + pinj
            + physics_variables.pohmmw
            - physics_variables.pradmw
        )

        # The following line is unphysical, but prevents -ve sqrt argument
        # Should be obsolete if constraint eqn 17 is turned on
        physics_variables.pdivt = max(0.001e0, physics_variables.pdivt)

        # if double null configuration share the power
        # over the upper and lower divertor, where physics_variables.ftar gives
        # the factor of power conducted to the lower divertor
        if physics_variables.idivrt == 2:
            physics_variables.pdivl = physics_variables.ftar * physics_variables.pdivt
            physics_variables.pdivu = (
                1.0e0 - physics_variables.ftar
            ) * physics_variables.pdivt
            physics_variables.pdivmax = max(
                physics_variables.pdivl, physics_variables.pdivu
            )

        # Resistive diffusion time = current penetration time ~ mu0.a^2/resistivity
        physics_variables.res_time = physics_functions_module.res_diff_time()

        # Power transported to the first wall by escaped alpha particles
        physics_variables.palpfwmw = physics_variables.palpmw * (
            1.0e0 - physics_variables.falpha
        )

        # Density limit
        physics_variables.dlimit, physics_variables.dnelimt = self.culdlm(
            physics_variables.bt,
            physics_variables.idensl,
            physics_variables.pdivt,
            physics_variables.plascur,
            divertor_variables.prn1,
            physics_variables.qstar,
            physics_variables.q95,
            physics_variables.rmajor,
            physics_variables.rminor,
            physics_variables.sarea,
            physics_variables.zeff,
        )

        # Calculate transport losses and energy confinement time using the
        # chosen scaling law
        (
            physics_variables.kappaa,
            physics_variables.ptrepv,
            physics_variables.ptripv,
            physics_variables.tauee,
            physics_variables.taueff,
            physics_variables.tauei,
            physics_variables.powerht,
        ) = self.pcond(
            physics_variables.afuel,
            physics_variables.palpmw,
            physics_variables.aspect,
            physics_variables.bt,
            physics_variables.dnitot,
            physics_variables.dene,
            physics_variables.dnla,
            physics_variables.eps,
            physics_variables.hfact,
            physics_variables.iinvqd,
            physics_variables.isc,
            physics_variables.ignite,
            physics_variables.kappa,
            physics_variables.kappa95,
            physics_variables.pchargemw,
            current_drive_variables.pinjmw,
            physics_variables.plascur,
            physics_variables.pcoreradpv,
            physics_variables.rmajor,
            physics_variables.rminor,
            physics_variables.te,
            physics_variables.ten,
            physics_variables.tin,
            physics_variables.q95,
            physics_variables.qstar,
            physics_variables.vol,
            physics_variables.xarea,
            physics_variables.zeff,
        )

        physics_variables.ptremw = physics_variables.ptrepv * physics_variables.vol
        physics_variables.ptrimw = physics_variables.ptripv * physics_variables.vol
        # Total transport power from scaling law (MW)
        # pscalingmw = physics_variables.ptremw + physics_variables.ptrimw #KE - why is this commented?

        # Calculate Volt-second requirements
        (
            physics_variables.phiint,
            physics_variables.rlp,
            physics_variables.vsbrn,
            physics_variables.vsind,
            physics_variables.vsres,
            physics_variables.vsstt,
        ) = vscalc(
            physics_variables.csawth,
            physics_variables.eps,
            physics_variables.facoh,
            physics_variables.gamma,
            physics_variables.kappa,
            physics_variables.rmajor,
            physics_variables.rplas,
            physics_variables.plascur,
            times_variables.t_fusion_ramp,
            times_variables.tburn,
            physics_variables.rli,
            constants.rmu0,
        )

        # Calculate auxiliary physics related information
        sbar = 1.0e0
        (
            physics_variables.burnup,
            physics_variables.dntau,
            physics_variables.figmer,
            physics_module.fusrat,
            physics_variables.qfuel,
            physics_variables.rndfuel,
            physics_variables.taup,
        ) = self.phyaux(
            physics_variables.aspect,
            physics_variables.dene,
            physics_variables.deni,
            physics_variables.fusionrate,
            physics_variables.alpharate,
            physics_variables.plascur,
            sbar,
            physics_variables.dnalp,
            physics_variables.taueff,
            physics_variables.vol,
        )

        # ptremw = physics_variables.ptrepv*physics_variables.vol
        # ptrimw = physics_variables.ptripv*physics_variables.vol
        # Total transport power from scaling law (MW)
        physics_variables.pscalingmw = (
            physics_variables.ptremw + physics_variables.ptrimw
        )

        # Calculate physics_variables.beta limit

        if physics_variables.iprofile == 1:
            # Relation between physics_variables.beta limit and plasma internal inductance
            # Hartmann and Zohm
            physics_variables.dnbeta = 4.0e0 * physics_variables.rli

        if physics_variables.iprofile == 2:
            # Original scaling law
            physics_variables.dnbeta = 2.7e0 * (
                1.0e0 + 5.0e0 * physics_variables.eps**3.5e0
            )

        if physics_variables.iprofile == 3 or physics_variables.iprofile == 5:
            # physics_variables.dnbeta found from physics_variables.aspect ratio scaling on p32 of Menard:
            # Menard, et al. "Fusion Nuclear Science Facilities
            # and Pilot Plants Based on the Spherical Tokamak."
            # Nucl. Fusion, 2016, 44.
            physics_variables.dnbeta = 3.12e0 + 3.5e0 * physics_variables.eps**1.7e0

        if physics_variables.iprofile == 6:
            # Method used for STEP plasma scoping
            # Tholerus et al. (2024), arXiv:2403.09460
            Fp = (physics_variables.ne0 * physics_variables.te0) / (
                physics_variables.dene * physics_variables.te
            )
            physics_variables.dnbeta = 3.7e0 + (
                (physics_variables.c_beta / Fp) * (12.5e0 - 3.5e0 * Fp)
            )

        # culblm returns the betalim for beta
        physics_variables.betalim = culblm(
            physics_variables.bt,
            physics_variables.dnbeta,
            physics_variables.plascur,
            physics_variables.rminor,
        )

        # MDK
        # Nominal mean photon wall load on entire first wall area including divertor and beam holes
        # Note that 'fwarea' excludes these, so they have been added back in.
        if physics_variables.iwalld == 1:
            physics_variables.photon_wall = (
                physics_variables.ffwal
                * physics_variables.pradmw
                / physics_variables.sarea
            )
        else:
            if physics_variables.idivrt == 2:
                # Double Null configuration in - including SoL radiation
                physics_variables.photon_wall = (
                    1.0e0 - fwbs_variables.fhcd - 2.0e0 * fwbs_variables.fdiv
                ) * physics_variables.pradmw / build_variables.fwarea + (
                    1.0e0 - fwbs_variables.fhcd - 2.0e0 * fwbs_variables.fdiv
                ) * physics_variables.rad_fraction_sol * physics_variables.pdivt / (
                    build_variables.fwarea
                )
            else:
                # Single null configuration - including SoL radaition
                physics_variables.photon_wall = (
                    1.0e0 - fwbs_variables.fhcd - fwbs_variables.fdiv
                ) * physics_variables.pradmw / build_variables.fwarea + (
                    1.0e0 - fwbs_variables.fhcd - fwbs_variables.fdiv
                ) * physics_variables.rad_fraction_sol * physics_variables.pdivt / build_variables.fwarea

        constraint_variables.peakradwallload = (
            physics_variables.photon_wall * constraint_variables.peakfactrad
        )

        # Calculate the target imbalances
        # find the total power into the targets
        physics_module.ptarmw = physics_variables.pdivt * (
            1.0e0 - physics_variables.rad_fraction_sol
        )
        # use physics_variables.ftar to find deltarsep
        # Parameters taken from double null machine
        # D. Brunner et al
        physics_module.lambdaio = 1.57e-3

        # Issue #1559 Infinities in physics_module.drsep when running single null in a double null machine
        # C W Ashe
        if physics_variables.ftar < 4.5e-5:
            physics_module.drsep = 1.5e-2
        elif physics_variables.ftar > (1.0e0 - 4.5e-5):
            physics_module.drsep = -1.5e-2
        else:
            physics_module.drsep = (
                -2.0e0 * 1.5e-3 * math.atanh(2.0e0 * (physics_variables.ftar - 0.5e0))
            )
        # Model Taken from D3-D paper for conventional divertor
        # Journal of Nuclear Materials
        # Volumes 290–293, March 2001, Pages 935-939
        # Find the innner and outer lower target imbalance

        physics_module.fio = 0.16e0 + (0.16e0 - 0.41e0) * (
            1.0e0
            - (
                2.0e0
                / (
                    1.0e0
                    + np.exp(-((physics_module.drsep / physics_module.lambdaio) ** 2))
                )
            )
        )
        if physics_variables.idivrt == 2:
            # Double Null configuration
            # Find all the power fractions accross the targets
            # Taken from D3-D conventional divertor design
            physics_module.fli = physics_variables.ftar * physics_module.fio
            physics_module.flo = physics_variables.ftar * (1.0e0 - physics_module.fio)
            physics_module.fui = (1.0e0 - physics_variables.ftar) * physics_module.fio
            physics_module.fuo = (1.0e0 - physics_variables.ftar) * (
                1.0e0 - physics_module.fio
            )
            # power into each target
            physics_module.plimw = physics_module.fli * physics_module.ptarmw
            physics_module.plomw = physics_module.flo * physics_module.ptarmw
            physics_module.puimw = physics_module.fui * physics_module.ptarmw
            physics_module.puomw = physics_module.fuo * physics_module.ptarmw
        else:
            # Single null configuration
            physics_module.fli = physics_module.fio
            physics_module.flo = 1.0e0 - physics_module.fio
            # power into each target
            physics_module.plimw = physics_module.fli * physics_module.ptarmw
            physics_module.plomw = physics_module.flo * physics_module.ptarmw

        # Calculate some derived quantities that may not have been defined earlier
        physics_module.total_loss_power = 1e6 * (
            physics_variables.falpha * physics_variables.palpmw
            + physics_variables.pchargemw
            + physics_variables.pohmmw
            + current_drive_variables.pinjmw
        )
        physics_module.rad_fraction_lcfs = (
            1.0e6 * physics_variables.pradmw / physics_module.total_loss_power
        )
        physics_variables.rad_fraction_total = (
            physics_module.rad_fraction_lcfs
            + (1.0e0 - physics_module.rad_fraction_lcfs)
            * physics_variables.rad_fraction_sol
        )
        physics_variables.pradsolmw = (
            physics_variables.rad_fraction_sol * physics_variables.pdivt
        )
        physics_module.total_plasma_internal_energy = (
            1.5e0
            * physics_variables.beta
            * physics_variables.btot
            * physics_variables.btot
            / (2.0e0 * constants.rmu0)
            * physics_variables.vol
        )

        physics_module.total_energy_conf_time = (
            physics_module.total_plasma_internal_energy
            / physics_module.total_loss_power
        )

        if any(numerics.icc == 78):
            po.write(
                self.outfile,
                (
                    f"reinke t and fz, physics = {physics_variables.tesep} , {reinke_variables.fzmin}"
                ),
            )
            # fsep = physics_variables.nesep / physics_variables.dene
            fgw = physics_variables.dlimit(7) / physics_variables.dene
            # calculate separatrix temperature, if Reinke criterion is used
            physics_variables.tesep = reinke_module.reinke_tsep(
                physics_variables.bt,
                constraint_variables.flhthresh,
                physics_variables.q95,
                physics_variables.rmajor,
                physics_variables.eps,
                fgw,
                physics_variables.kappa,
                reinke_variables.lhat,
            )

            if reinke_variables.fzmin >= 1.0e0:
                error_handling.report_error(217)

            po.write(
                self.outfile,
                (
                    f" 'fzactual, frac, reinke_variables.impvardiv = {reinke_variables.fzactual},"
                    f" {impurity_radiation_module.impurity_arr_frac(reinke_variables.impvardiv)},"
                    f" {reinke_variables.impvardiv}"
                ),
            )

    @staticmethod
    def culdlm(
        bt, idensl, pdivt, plascur, prn1, qcyl, q95, rmajor, rminor, sarea, zeff
    ):
        """Density limit calculation
        author: P J Knight, CCFE, Culham Science Centre
        bt       : input real :  toroidal field on axis (T)
        idensl   : input/output integer : switch denoting which formula to enforce
        pdivt    : input real :  power flowing to the edge plasma via
        charged particles (MW)
        plascur  : input real :  plasma current (A)
        prn1     : input real :  edge density / average plasma density
        qcyl     : input real :  equivalent cylindrical safety factor (qstar)
        q95      : input real :  safety factor at 95% surface
        rmajor   : input real :  plasma major radius (m)
        rminor   : input real :  plasma minor radius (m)
        sarea    : input real :  plasma surface area (m**2)
        zeff     : input real :  plasma effective charge
        dlimit(7): output real array : average plasma density limit using
        seven different models (m**-3)
        dnelimt  : output real : enforced average plasma density limit (m**-3)
        This routine calculates several different formulae for the
        density limit, and enforces the one chosen by the user.
        AEA FUS 172: Physics Assessment for the European Reactor Study
        """
        if idensl < 1 or idensl > 7:
            error_handling.idiags[0] = idensl
            error_handling.report_error(79)

        dlimit = np.empty((7,))

        # Power per unit area crossing the plasma edge
        # (excludes radiation and neutrons)

        qperp = pdivt / sarea

        # Old ASDEX density limit formula
        # This applies to the density at the plasma edge, so must be scaled
        # to give the density limit applying to the average plasma density.

        dlimit[0] = (
            1.54e20 * qperp**0.43 * bt**0.31 / (q95 * rmajor) ** 0.45
        ) / prn1

        # Borrass density limit model for ITER (I)
        # This applies to the density at the plasma edge, so must be scaled
        # to give the density limit applying to the average plasma density.
        # Borrass et al, ITER-TN-PH-9-6 (1989)

        dlimit[1] = (
            1.8e20 * qperp**0.53 * bt**0.31 / (q95 * rmajor) ** 0.22
        ) / prn1

        # Borrass density limit model for ITER (II)
        # This applies to the density at the plasma edge, so must be scaled
        # to give the density limit applying to the average plasma density.
        # This formula is (almost) identical to that in the original routine
        # denlim (now deleted).

        dlimit[2] = (
            0.5e20 * qperp**0.57 * bt**0.31 / (q95 * rmajor) ** 0.09
        ) / prn1

        # JET edge radiation density limit model
        # This applies to the density at the plasma edge, so must be scaled
        # to give the density limit applying to the average plasma density.
        # qcyl=qstar here, but literature is not clear.

        denom = (zeff - 1.0) * (1.0 - 4.0 / (3.0 * qcyl))
        if denom <= 0.0:
            if idensl == 4:
                error_handling.fdiags[0] = denom
                error_handling.fdiags[1] = qcyl
                error_handling.report_error(80)
                idensl = 5

            dlimit[3] = 0.0
        else:
            dlimit[3] = (1.0e20 * np.sqrt(pdivt / denom)) / prn1

        # JET simplified density limit model
        # This applies to the density at the plasma edge, so must be scaled
        # to give the density limit applying to the average plasma density.

        dlimit[4] = (0.237e20 * bt * np.sqrt(pdivt) / rmajor) / prn1

        # Hugill-Murakami M.q limit
        # qcyl=qstar here, which is okay according to the literature

        dlimit[5] = 3.0e20 * bt / (rmajor * qcyl)

        # Greenwald limit

        dlimit[6] = 1.0e14 * plascur / (np.pi * rminor * rminor)

        # Enforce the chosen density limit

        return dlimit, dlimit[idensl - 1]

    @staticmethod
    def plasma_composition():
        """Calculates various plasma component fractional makeups
        author: P J Knight, CCFE, Culham Science Centre

        This subroutine determines the various plasma component
        fractional makeups. It is the replacement for the original
        It is the replacement for the original routine <CODE>betcom</CODE>,
        and is used in conjunction with the new impurity radiation model
        """
        # Alpha ash portion
        physics_variables.dnalp = physics_variables.dene * physics_variables.ralpne

        # Protons
        # This calculation will be wrong on the first call as the particle
        # production rates are evaluated later in the calling sequence
        # Issue #557 Allow protium impurity to be specified: 'protium'
        # This will override the calculated value which is a minimum.
        if physics_variables.alpharate < 1.0e-6:  # not calculated yet...
            physics_variables.dnprot = max(
                physics_variables.protium * physics_variables.dene,
                physics_variables.dnalp * (physics_variables.fhe3 + 1.0e-3),
            )  # rough estimate
        else:
            physics_variables.dnprot = max(
                physics_variables.protium * physics_variables.dene,
                physics_variables.dnalp
                * physics_variables.protonrate
                / physics_variables.alpharate,
            )

        # Beam hot ion component
        # If ignited, prevent beam fusion effects
        if physics_variables.ignite == 0:
            physics_variables.dnbeam = physics_variables.dene * physics_variables.rnbeam
        else:
            physics_variables.dnbeam = 0.0

        # Sum of Zi.ni for all impurity ions (those with charge > helium)
        znimp = 0.0
        for imp in range(impurity_radiation_module.nimp):
            if impurity_radiation_module.impurity_arr_z[imp] > 2:
                znimp += impurity_radiation_module.zav_of_te(
                    imp + 1, physics_variables.te
                ) * (
                    impurity_radiation_module.impurity_arr_frac[imp]
                    * physics_variables.dene
                )

        # Fuel portion - conserve charge neutrality
        # znfuel is the sum of Zi.ni for the three fuel ions
        znfuel = (
            physics_variables.dene
            - 2.0 * physics_variables.dnalp
            - physics_variables.dnprot
            - physics_variables.dnbeam
            - znimp
        )

        # Fuel ion density, deni
        # For D-T-He3 mix, deni = nD + nT + nHe3, while znfuel = nD + nT + 2*nHe3
        # So deni = znfuel - nHe3 = znfuel - fhe3*deni
        physics_variables.deni = znfuel / (1.0 + physics_variables.fhe3)

        # Set hydrogen and helium impurity fractions for
        # radiation calculations
        impurity_radiation_module.impurity_arr_frac[
            impurity_radiation_module.element2index("H_") - 1
        ] = (
            physics_variables.dnprot
            + (physics_variables.fdeut + physics_variables.ftrit)
            * physics_variables.deni
            + physics_variables.dnbeam
        ) / physics_variables.dene

        impurity_radiation_module.impurity_arr_frac[
            impurity_radiation_module.element2index("He") - 1
        ] = (
            physics_variables.fhe3 * physics_variables.deni / physics_variables.dene
            + physics_variables.ralpne
        )

        # Total impurity density
        physics_variables.dnz = 0.0
        for imp in range(impurity_radiation_module.nimp):
            if impurity_radiation_module.impurity_arr_z[imp] > 2:
                physics_variables.dnz += (
                    impurity_radiation_module.impurity_arr_frac[imp]
                    * physics_variables.dene
                )

        # Total ion density
        physics_variables.dnitot = (
            physics_variables.deni
            + physics_variables.dnalp
            + physics_variables.dnprot
            + physics_variables.dnbeam
            + physics_variables.dnz
        )

        # Set some (obsolescent) impurity fraction variables
        # for the benefit of other routines
        physics_variables.rncne = impurity_radiation_module.impurity_arr_frac[
            impurity_radiation_module.element2index("C_")
        ]
        physics_variables.rnone = impurity_radiation_module.impurity_arr_frac[
            impurity_radiation_module.element2index("O_")
        ]
        physics_variables.rnfene = (
            impurity_radiation_module.impurity_arr_frac[
                impurity_radiation_module.element2index("Fe")
            ]
            + impurity_radiation_module.impurity_arr_frac[
                impurity_radiation_module.element2index("Ar")
            ]
        )

        # Effective charge
        # Calculation should be sum(ni.Zi^2) / sum(ni.Zi),
        # but ne = sum(ni.Zi) through quasineutrality
        physics_variables.zeff = 0.0
        for imp in range(impurity_radiation_module.nimp):
            physics_variables.zeff += (
                impurity_radiation_module.impurity_arr_frac[imp]
                * impurity_radiation_module.zav_of_te(imp + 1, physics_variables.te)
                ** 2
            )

        # Define coulomb logarithm
        # (collisions: ion-electron, electron-electron)
        physics_variables.dlamee = (
            31.0
            - (np.log(physics_variables.dene) / 2.0)
            + np.log(physics_variables.te * 1000.0)
        )
        physics_variables.dlamie = (
            31.3
            - (np.log(physics_variables.dene) / 2.0)
            + np.log(physics_variables.te * 1000.0)
        )

        # Fraction of alpha energy to ions and electrons
        # From Max Fenstermacher
        # (used with electron and ion power balance equations only)
        # No consideration of pchargepv here...

        # pcoef now calculated in plasma_profiles, after the very first
        # call of plasma_composition; use old parabolic profile estimate
        # in this case
        if physics_module.first_call == 1:
            pc = (
                (1.0 + physics_variables.alphan)
                * (1.0 + physics_variables.alphat)
                / (1.0 + physics_variables.alphan + physics_variables.alphat)
            )
            physics_module.first_call = 0
        else:
            pc = physics_variables.pcoef

        physics_variables.falpe = 0.88155 * np.exp(-physics_variables.te * pc / 67.4036)
        physics_variables.falpi = 1.0 - physics_variables.falpe

        # Average atomic masses
        physics_variables.afuel = (
            2.0 * physics_variables.fdeut
            + 3.0 * physics_variables.ftrit
            + 3.0 * physics_variables.fhe3
        )
        physics_variables.abeam = (
            2.0 * (1.0 - current_drive_variables.ftritbm)
            + 3.0 * current_drive_variables.ftritbm
        )

        # Density weighted mass
        physics_variables.aion = (
            physics_variables.afuel * physics_variables.deni
            + 4.0 * physics_variables.dnalp
            + physics_variables.dnprot
            + physics_variables.abeam * physics_variables.dnbeam
        )
        for imp in range(impurity_radiation_module.nimp):
            if impurity_radiation_module.impurity_arr_z[imp] > 2:
                physics_variables.aion += (
                    physics_variables.dene
                    * impurity_radiation_module.impurity_arr_frac[imp]
                    * impurity_radiation_module.impurity_arr_amass[imp]
                )

        physics_variables.aion = physics_variables.aion / physics_variables.dnitot

        # Mass weighted plasma effective charge
        physics_variables.zeffai = (
            physics_variables.fdeut * physics_variables.deni / 2.0
            + physics_variables.ftrit * physics_variables.deni / 3.0
            + 4.0 * physics_variables.fhe3 * physics_variables.deni / 3.0
            + physics_variables.dnalp
            + physics_variables.dnprot
            + (1.0 - current_drive_variables.ftritbm) * physics_variables.dnbeam / 2.0
            + current_drive_variables.ftritbm * physics_variables.dnbeam / 3.0
        ) / physics_variables.dene
        for imp in range(impurity_radiation_module.nimp):
            if impurity_radiation_module.impurity_arr_z[imp] > 2:
                physics_variables.zeffai += (
                    impurity_radiation_module.impurity_arr_frac[imp]
                    * impurity_radiation_module.zav_of_te(imp + 1, physics_variables.te)
                    ** 2
                    / impurity_radiation_module.impurity_arr_amass[imp]
                )

    @staticmethod
    def phyaux(
        aspect,
        dene,
        deni,
        fusionrate,
        alpharate,
        plascur,
        sbar,
        dnalp,
        taueff,
        vol,
    ):
        """Auxiliary physics quantities
        author: P J Knight, CCFE, Culham Science Centre
        aspect : input real :  plasma aspect ratio
        dene   : input real :  electron density (/m3)
        deni   : input real :  fuel ion density (/m3)
        dnalp  : input real :  alpha ash density (/m3)
        fusionrate : input real :  fusion reaction rate (/m3/s)
        alpharate  : input real :  alpha particle production rate (/m3/s)
        plascur: input real :  plasma current (A)
        sbar   : input real :  exponent for aspect ratio (normally 1)
        taueff : input real :  global energy confinement time (s)
        vol    : input real :  plasma volume (m3)
        burnup : output real : fractional plasma burnup
        dntau  : output real : plasma average n-tau (s/m3)
        figmer : output real : physics figure of merit
        fusrat : output real : number of fusion reactions per second
        qfuel  : output real : fuelling rate for D-T (nucleus-pairs/sec)
        rndfuel: output real : fuel burnup rate (reactions/s)
        taup   : output real : (alpha) particle confinement time (s)
        This subroutine calculates extra physics related items
        needed by other parts of the code
        """

        figmer = 1e-6 * plascur * aspect**sbar

        dntau = taueff * dene

        # Fusion reactions per second

        fusrat = fusionrate * vol

        # Alpha particle confinement time (s)
        # Number of alphas / alpha production rate

        if alpharate != 0.0:
            taup = dnalp / alpharate
        else:  # only likely if DD is only active fusion reaction
            taup = 0.0

        # Fractional burnup

        # (Consider detailed model in: G. L. Jackson, V. S. Chan, R. D. Stambaugh,
        # Fusion Science and Technology, vol.64, no.1, July 2013, pp.8-12)

        # The ratio of ash to fuel particle confinement times is given by
        # tauratio
        # Possible logic...
        # burnup = fuel ion-pairs burned/m3 / initial fuel ion-pairs/m3;
        # fuel ion-pairs burned/m3 = alpha particles/m3 (for both D-T and D-He3 reactions)
        # initial fuel ion-pairs/m3 = burnt fuel ion-pairs/m3 + unburnt fuel-ion pairs/m3
        # Remember that unburnt fuel-ion pairs/m3 = 0.5 * unburnt fuel-ions/m3
        if physics_variables.burnup_in <= 1.0e-9:
            burnup = dnalp / (dnalp + 0.5 * deni) / physics_variables.tauratio
        else:
            burnup = physics_variables.burnup_in
        # Fuel burnup rate (reactions/second) (previously Amps)

        rndfuel = fusrat

        # Required fuelling rate (fuel ion pairs/second) (previously Amps)

        qfuel = rndfuel / burnup

        return burnup, dntau, figmer, fusrat, qfuel, rndfuel, taup

    @staticmethod
    def pohm(facoh, kappa95, plascur, rmajor, rminor, ten, vol, zeff):
        # Density weighted electron temperature in 10 keV units

        t10 = ten / 10.0

        # Plasma resistance, from loop voltage calculation in IPDG89

        rplas = (
            physics_variables.plasma_res_factor
            * 2.15e-9
            * zeff
            * rmajor
            / (kappa95 * rminor**2 * t10**1.5)
        )

        # Neo-classical resistivity enhancement factor
        # Taken from  N. A. Uckan et al, Fusion Technology 13 (1988) p.411.
        # The expression is valid for aspect ratios in the range 2.5--4.

        rpfac = 4.3 - 0.6 * rmajor / rminor
        rplas = rplas * rpfac

        # Check to see if plasma resistance is negative
        # (possible if aspect ratio is too high)

        if rplas <= 0.0:
            error_handling.fdiags[0] = rplas
            error_handling.fdiags[1] = physics_variables.aspect
            error_handling.report_error(83)

        # Ohmic heating power per unit volume
        # Corrected from: pohmpv = (facoh*plascur)**2 * ...

        pohmpv = facoh * plascur**2 * rplas * 1.0e-6 / vol

        # Total ohmic heating power

        pohmmw = pohmpv * vol

        return pohmpv, pohmmw, rpfac, rplas

    @staticmethod
    def culcur(
        alphaj,
        alphap,
        bt,
        eps,
        icurr,
        iprofile,
        kappa,
        kappa95,
        p0,
        pperim,
        q0,
        q95,
        rli,
        rmajor,
        rminor,
        sf,
        triang,
        triang95,
    ):
        """Routine to calculate the plasma current
        author: P J Knight, CCFE, Culham Science Centre
        alphaj   : input/output real : current profile index
        alphap   : input real :  pressure profile index
        bt       : input real :  toroidal field on axis (T)
        eps      : input real :  inverse aspect ratio
        icurr    : input integer : current scaling model to use
        1 = Peng analytic fit
        2 = Peng divertor scaling (TART)
        3 = simple ITER scaling
        4 = revised ITER scaling
        5 = Todd empirical scaling I
        6 = Todd empirical scaling II
        7 = Connor-Hastie model
        8 = Sauter scaling (allowing negative triangularity) Issue #392
            'Geometric formulas for system codes including the effect of negative triangularity'
        iprofile : input integer : switch for current profile consistency
        0 use input values for alphaj, rli, dnbeta
        1 make these consistent with input q, q_0 values (recommend `icurr=4` with this option)
        2 use input values for alphaj, rli. Scale dnbeta with aspect ratio (original scaling)
        3 use input values for alphaj, rli. Scale dnbeta with aspect ratio (Menard scaling)
        4 use input values for alphaj, dnbeta. Set rli from elongation (Menard scaling)
        5 use input value for alphaj.  Set rli and dnbeta from Menard scaling
        kappa    : input real :  plasma elongation
        kappa95  : input real :  plasma elongation at 95% surface
        p0       : input real :  central plasma pressure (Pa)
        pperim   : input real :  plasma perimeter length (m)
        q0       : input real :  plasma safety factor on axis
        q95      : input real :  plasma safety factor at 95% flux (= q-bar for icurr=2)
        rli      : input/output real : plasma normalised internal inductance
        rmajor   : input real :  major radius (m)
        rminor   : input real :  minor radius (m)
        sf       : input real :  shape factor for icurr=1 (=A/pi in documentation)
        triang   : input real :  plasma triangularity
        triang95 : input real :  plasma triangularity at 95% surface
        bp       : output real : poloidal field (T)
        qstar    : output real : equivalent cylindrical safety factor (shaped)
        plascur  : output real : plasma current (A)
        This routine calculates the plasma current based on the edge
        safety factor q95. It will also make the current profile parameters
        consistent with the q-profile if required.
        J D Galambos, STAR Code : Spherical Tokamak Analysis and Reactor Code,
        unpublished internal Oak Ridge document
        Y.-K. M. Peng, J. Galambos and P.C. Shipe, 1992,
        Fusion Technology, 21, 1729
        ITER Physics Design Guidelines: 1989 [IPDG89], N. A. Uckan et al,
        ITER Documentation Series No.10, IAEA/ITER/DS/10, IAEA, Vienna, 1990
        M. Kovari et al, 2014, "PROCESS": A systems code for fusion power plants -
        Part 1: Physics https://www.sciencedirect.com/science/article/pii/S0920379614005961
        H. Zohm et al, 2013, On the Physics Guidelines for a Tokamak DEMO
        https://iopscience.iop.org/article/10.1088/0029-5515/53/7/073019
        T. Hartmann, 2013, Development of a modular systems code to analyse the
        implications of physics assumptions on the design of a demonstration fusion power plant
        https://inis.iaea.org/search/search.aspx?orig_q=RN:45031642
        Sauter, Geometric formulas for systems codes..., FED 2016
        """
        # Aspect ratio

        asp = 1.0 / eps

        # Calculate the function Fq that scales the edge q from the
        # circular cross-section cylindrical case

        # Only the Sauter scaling (icurr=8) is suitable for negative triangularity:

        if icurr != 8 and triang < 0.0:
            raise ValueError(
                f"Triangularity is negative without icurr = 8: {triang=}, {icurr=}"
            )

        if icurr == 1:  # Peng analytical fit
            fq = (1.22 - 0.68 * eps) / ((1.0 - eps * eps) ** 2) * sf**2
        elif icurr == 2:  # Peng scaling for double null divertor; TARTs [STAR Code]
            plascur = 1.0e6 * plasc(q95, asp, eps, rminor, bt, kappa, triang)
        elif icurr == 3:  # Simple ITER scaling (simply the cylindrical case)
            fq = 1.0
        elif icurr == 4:  # ITER formula (IPDG89)
            fq = (
                0.5
                * (1.17 - 0.65 * eps)
                / ((1.0 - eps * eps) ** 2)
                * (
                    1.0
                    + kappa95**2 * (1.0 + 2.0 * triang95**2 - 1.2 * triang95**3)
                )
            )
        elif icurr in [5, 6]:  # Todd empirical scalings
            fq = (
                (1.0 + 2.0 * eps * eps)
                * 0.5
                * (1.0 + kappa95**2)
                * (
                    1.24
                    - 0.54 * kappa95
                    + 0.3 * (kappa95**2 + triang95**2)
                    + 0.125 * triang95
                )
            )

            fq *= 1 if icurr == 7 else (1.0 + (abs(kappa95 - 1.2)) ** 3)
        elif icurr == 7:  # Connor-Hastie asymptotically-correct expression
            # N.B. If iprofile=1, alphaj will be wrong during the first call (only)
            fq = conhas(alphaj, alphap, bt, triang95, eps, kappa95, p0, constants.rmu0)
        elif (
            icurr == 8
        ):  # Sauter scaling allowing negative triangularity [FED May 2016]
            # Assumes zero squareness, note takes kappa, delta at separatrix not _95

            w07 = 1.0  # zero squareness - can be modified later if required

            fq = (
                (4.1e6 / 5.0e6)
                * (1.0 + 1.2 * (kappa - 1.0) + 0.56 * (kappa - 1.0) ** 2)
                * (1.0 + 0.09 * triang + 0.16 * triang**2)
                * (1.0 + 0.45 * triang * eps)
                / (1.0 - 0.74 * eps)
                * (1.0 + 0.55 * (w07 - 1.0))
            )
        elif icurr == 9:
            fq = 0.538 * (1.0 + 2.440 * eps**2.736) * kappa**2.154 * triang**0.060
        else:
            raise ValueError(f"Invalid value {icurr=}")

        if icurr != 2:
            plascur = 5.0e6 * rminor**2 / (rmajor * q95) * fq * bt
        # == 2 case covered above

        qstar = (
            5.0e6
            * rminor**2
            / (rmajor * plascur / bt)
            * 0.5
            * (1.0 + kappa95**2 * (1.0 + 2.0 * triang95**2 - 1.2 * triang95**3))
        )

        physics_variables.normalised_total_beta = (
            1.0e8 * physics_variables.beta * rminor * bt / plascur
        )
        bp = bpol(
            icurr, plascur, q95, asp, eps, bt, kappa, triang, pperim, constants.rmu0
        )

        if iprofile == 1:
            # Ensure current profile consistency, if required
            # This is as described in Hartmann and Zohm only if icurr = 4 as well...

            # Tokamaks 4th Edition, Wesson, page 116
            alphaj = qstar / q0 - 1.0
            rli = np.log(1.65 + 0.89 * alphaj)

        if iprofile in [4, 5, 6]:
            # Spherical Tokamak relation for internal inductance
            # Menard et al. (2016), Nuclear Fusion, 56, 106023
            rli = 3.4 - kappa

        return alphaj, rli, bp, qstar, plascur

    def outtim(self):
        po.oheadr(self.outfile, "Times")

        po.ovarrf(
            self.outfile,
            "Initial charge time for CS from zero current (s)",
            "(tramp)",
            times_variables.tramp,
        )
        po.ovarrf(
            self.outfile,
            "Plasma current ramp-up time (s)",
            "(tohs)",
            times_variables.tohs,
        )
        po.ovarrf(
            self.outfile,
            "Heating time (s)",
            "(t_fusion_ramp)",
            times_variables.t_fusion_ramp,
        )
        po.ovarre(
            self.outfile, "Burn time (s)", "(tburn)", times_variables.tburn, "OP "
        )
        po.ovarrf(
            self.outfile,
            "Reset time to zero current for CS (s)",
            "(tqnch)",
            times_variables.tqnch,
        )
        po.ovarrf(
            self.outfile, "Time between pulses (s)", "(tdwell)", times_variables.tdwell
        )
        po.oblnkl(self.outfile)
        po.ovarre(
            self.outfile,
            "Total plant cycle time (s)",
            "(tcycle)",
            times_variables.tcycle,
            "OP ",
        )

    def outplas(self):
        """Subroutine to output the plasma physics information
        author: P J Knight, CCFE, Culham Science Centre
        self.outfile : input integer : Fortran output unit identifier
        This routine writes the plasma physics information
        to a file, in a tidy format.
        """

        # ###############################################
        # Dimensionless plasma parameters. See reference below.
        physics_module.nu_star = (
            1
            / constants.rmu0
            * (15.0e0 * constants.echarge**4 * physics_variables.dlamie)
            / (4.0e0 * np.pi**1.5e0 * constants.epsilon0**2)
            * physics_variables.vol**2
            * physics_variables.rmajor**2
            * physics_variables.bt
            * np.sqrt(physics_variables.eps)
            * physics_variables.dnla**3
            * physics_variables.kappa
            / (
                physics_module.total_plasma_internal_energy**2
                * physics_variables.plascur
            )
        )

        physics_module.rho_star = np.sqrt(
            2.0e0
            * constants.mproton
            * physics_variables.aion
            * physics_module.total_plasma_internal_energy
            / (3.0e0 * physics_variables.vol * physics_variables.dnla)
        ) / (
            constants.echarge
            * physics_variables.bt
            * physics_variables.eps
            * physics_variables.rmajor
        )

        physics_module.beta_mcdonald = (
            4.0e0
            / 3.0e0
            * constants.rmu0
            * physics_module.total_plasma_internal_energy
            / (physics_variables.vol * physics_variables.bt**2)
        )

        po.oheadr(self.outfile, "Plasma")

        if stellarator_variables.istell == 0:
            if physics_variables.idivrt == 0:
                po.ocmmnt(self.outfile, "Plasma configuration = limiter")
            elif physics_variables.idivrt == 1:
                po.ocmmnt(self.outfile, "Plasma configuration = single null divertor")
            elif physics_variables.idivrt == 2:
                po.ocmmnt(self.outfile, "Plasma configuration = double null divertor")
            else:
                error_handling.idiags[0] = physics_variables.idivrt
                po.report_error(85)
        else:
            po.ocmmnt(self.outfile, "Plasma configuration = stellarator")

        if stellarator_variables.istell == 0:
            if physics_variables.itart == 0:
                physics_module.itart_r = physics_variables.itart
                po.ovarrf(
                    self.outfile,
                    "Tokamak aspect ratio = Conventional, itart = 0",
                    "(itart)",
                    physics_module.itart_r,
                )
            elif physics_variables.itart == 1:
                physics_module.itart_r = physics_variables.itart
                po.ovarrf(
                    self.outfile,
                    "Tokamak aspect ratio = Spherical, itart = 1",
                    "(itart)",
                    physics_module.itart_r,
                )

        po.osubhd(self.outfile, "Plasma Geometry :")
        po.ovarrf(
            self.outfile, "Major radius (m)", "(rmajor)", physics_variables.rmajor
        )
        po.ovarrf(
            self.outfile,
            "Minor radius (m)",
            "(rminor)",
            physics_variables.rminor,
            "OP ",
        )
        po.ovarrf(self.outfile, "Aspect ratio", "(aspect)", physics_variables.aspect)

        if stellarator_variables.istell == 0:
            if physics_variables.ishape in [0, 6, 8]:
                po.ovarrf(
                    self.outfile,
                    "Elongation, X-point (input value used)",
                    "(kappa)",
                    physics_variables.kappa,
                    "IP ",
                )
            elif physics_variables.ishape == 1:
                po.ovarrf(
                    self.outfile,
                    "Elongation, X-point (TART scaling)",
                    "(kappa)",
                    physics_variables.kappa,
                    "OP ",
                )
            elif physics_variables.ishape in [2, 3]:
                po.ovarrf(
                    self.outfile,
                    "Elongation, X-point (Zohm scaling)",
                    "(kappa)",
                    physics_variables.kappa,
                    "OP ",
                )
                po.ovarrf(
                    self.outfile,
                    "Zohm scaling adjustment factor",
                    "(fkzohm)",
                    physics_variables.fkzohm,
                )
            elif physics_variables.ishape in [4, 5, 7]:
                po.ovarrf(
                    self.outfile,
                    "Elongation, X-point (calculated from kappa95)",
                    "(kappa)",
                    physics_variables.kappa,
                    "OP ",
                )
            elif physics_variables.ishape == 9:
                po.ovarrf(
                    self.outfile,
                    "Elongation, X-point (calculated from aspect ratio and li(3))",
                    "(kappa)",
                    physics_variables.kappa,
                    "OP ",
                )
            elif physics_variables.ishape == 10:
                po.ovarrf(
                    self.outfile,
                    "Elongation, X-point (calculated from aspect ratio and stability margin)",
                    "(kappa)",
                    physics_variables.kappa,
                    "OP ",
                )
            elif physics_variables.ishape == 11:
                po.ovarrf(
                    self.outfile,
                    "Elongation, X-point (calculated from aspect ratio via Menard 2016)",
                    "(kappa)",
                    physics_variables.kappa,
                    "OP ",
                )
            else:
                error_handling.idiags[0] = physics_variables.ishape
                po.report_error(86)

            if physics_variables.ishape in [4, 5, 7]:
                po.ovarrf(
                    self.outfile,
                    "Elongation, 95% surface (input value used)",
                    "(kappa95)",
                    physics_variables.kappa95,
                    "IP ",
                )
            else:
                po.ovarrf(
                    self.outfile,
                    "Elongation, 95% surface (calculated from kappa)",
                    "(kappa95)",
                    physics_variables.kappa95,
                    "OP ",
                )

            po.ovarrf(
                self.outfile,
                "Elongation, area ratio calc.",
                "(kappaa)",
                physics_variables.kappaa,
                "OP ",
            )

            if physics_variables.ishape in [0, 2, 6, 8, 9, 10, 11]:
                po.ovarrf(
                    self.outfile,
                    "Triangularity, X-point (input value used)",
                    "(triang)",
                    physics_variables.triang,
                    "IP ",
                )
            elif physics_variables.ishape == 1:
                po.ovarrf(
                    self.outfile,
                    "Triangularity, X-point (TART scaling)",
                    "(triang)",
                    physics_variables.triang,
                    "OP ",
                )
            else:
                po.ovarrf(
                    self.outfile,
                    "Triangularity, X-point (calculated from triang95)",
                    "(triang)",
                    physics_variables.triang,
                    "OP ",
                )

            if physics_variables.ishape in [3, 4, 5, 7]:
                po.ovarrf(
                    self.outfile,
                    "Triangularity, 95% surface (input value used)",
                    "(triang95)",
                    physics_variables.triang95,
                    "IP ",
                )
            else:
                po.ovarrf(
                    self.outfile,
                    "Triangularity, 95% surface (calculated from triang)",
                    "(triang95)",
                    physics_variables.triang95,
                    "OP ",
                )

            po.ovarrf(
                self.outfile,
                "Plasma poloidal perimeter (m)",
                "(pperim)",
                physics_variables.pperim,
                "OP ",
            )

            po.ovarre(
                self.outfile,
                "Plasma cross-sectional area (m2)",
                "(xarea)",
                physics_variables.xarea,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Plasma surface area (m2)",
                "(sarea)",
                physics_variables.sarea,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Plasma volume (m3)",
                "(vol)",
                physics_variables.vol,
                "OP ",
            )

            po.osubhd(self.outfile, "Current and Field :")

            if stellarator_variables.istell == 0:
                if physics_variables.iprofile == 1:
                    po.ocmmnt(
                        self.outfile,
                        "Consistency between q0,q,alphaj,rli,dnbeta is enforced",
                    )
                else:
                    po.ocmmnt(
                        self.outfile,
                        "Consistency between q0,q,alphaj,rli,dnbeta is not enforced",
                    )

                po.oblnkl(self.outfile)
                po.ovarin(
                    self.outfile,
                    "Plasma current scaling law used",
                    "(icurr)",
                    physics_variables.icurr,
                )

                po.ovarrf(
                    self.outfile,
                    "Plasma current (MA)",
                    "(plascur/1D6)",
                    physics_variables.plascur / 1.0e6,
                    "OP ",
                )

                if physics_variables.iprofile == 1:
                    po.ovarrf(
                        self.outfile,
                        "Current density profile factor",
                        "(alphaj)",
                        physics_variables.alphaj,
                        "OP ",
                    )
                else:
                    po.ovarrf(
                        self.outfile,
                        "Current density profile factor",
                        "(alphaj)",
                        physics_variables.alphaj,
                    )

                po.ovarrf(
                    self.outfile,
                    "Plasma internal inductance, li",
                    "(rli)",
                    physics_variables.rli,
                    "OP ",
                )
                po.ovarrf(
                    self.outfile,
                    "Vertical field at plasma (T)",
                    "(bvert)",
                    physics_variables.bvert,
                    "OP ",
                )

            po.ovarrf(
                self.outfile,
                "Vacuum toroidal field at R (T)",
                "(bt)",
                physics_variables.bt,
            )
            po.ovarrf(
                self.outfile,
                "Average poloidal field (T)",
                "(bp)",
                physics_variables.bp,
                "OP ",
            )

            po.ovarrf(
                self.outfile,
                "Total field (sqrt(bp^2 + bt^2)) (T)",
                "(btot)",
                physics_variables.btot,
                "OP ",
            )

        if stellarator_variables.istell == 0:
            po.ovarrf(
                self.outfile, "Safety factor on axis", "(q0)", physics_variables.q0
            )

            if physics_variables.icurr == 2:
                po.ovarrf(
                    self.outfile, "Mean edge safety factor", "(q)", physics_variables.q
                )

            po.ovarrf(
                self.outfile,
                "Safety factor at 95% flux surface",
                "(q95)",
                physics_variables.q95,
            )

            po.ovarrf(
                self.outfile,
                "Cylindrical safety factor (qcyl)",
                "(qstar)",
                physics_variables.qstar,
                "OP ",
            )

            if physics_variables.ishape == 1:
                po.ovarrf(
                    self.outfile,
                    "Lower limit for edge safety factor q",
                    "(qlim)",
                    physics_variables.qlim,
                    "OP ",
                )

        else:
            po.ovarrf(
                self.outfile,
                "Rotational transform",
                "(iotabar)",
                stellarator_variables.iotabar,
            )

        po.osubhd(self.outfile, "Beta Information :")

        betath = (
            physics_variables.beta - physics_variables.betaft - physics_variables.betanb
        )
        gammaft = (physics_variables.betaft + physics_variables.betanb) / betath

        po.ovarre(self.outfile, "Total plasma beta", "(beta)", physics_variables.beta)
        po.ovarre(
            self.outfile,
            "Total poloidal beta",
            "(betap)",
            physics_variables.betap,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Total toroidal beta",
            " ",
            physics_variables.beta
            * (physics_variables.btot / physics_variables.bt) ** 2,
            "OP ",
        )
        po.ovarre(
            self.outfile, "Fast alpha beta", "(betaft)", physics_variables.betaft, "OP "
        )
        po.ovarre(
            self.outfile, "Beam ion beta", "(betanb)", physics_variables.betanb, "OP "
        )
        po.ovarre(
            self.outfile,
            "(Fast alpha + beam physics_variables.beta)/(thermal physics_variables.beta)",
            "(gammaft)",
            gammaft,
            "OP ",
        )

        po.ovarre(self.outfile, "Thermal beta", " ", betath, "OP ")
        po.ovarre(
            self.outfile,
            "Thermal poloidal beta",
            " ",
            betath * (physics_variables.btot / physics_variables.bp) ** 2,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Thermal toroidal physics_variables.beta (= beta-exp)",
            " ",
            betath * (physics_variables.btot / physics_variables.bt) ** 2,
            "OP ",
        )

        po.ovarrf(
            self.outfile,
            "2nd stability physics_variables.beta : beta_p / (R/a)",
            "(eps*betap)",
            physics_variables.eps * physics_variables.betap,
            "OP ",
        )
        po.ovarrf(
            self.outfile,
            "2nd stability physics_variables.beta upper limit",
            "(epbetmax)",
            physics_variables.epbetmax,
        )

        if stellarator_variables.istell == 0:
            if physics_variables.iprofile == 1:
                po.ovarrf(
                    self.outfile,
                    "Beta g coefficient",
                    "(dnbeta)",
                    physics_variables.dnbeta,
                    "OP ",
                )
            else:
                po.ovarrf(
                    self.outfile,
                    "Beta g coefficient",
                    "(dnbeta)",
                    physics_variables.dnbeta,
                )

            po.ovarrf(
                self.outfile,
                "Normalised thermal beta",
                " ",
                1.0e8
                * betath
                * physics_variables.rminor
                * physics_variables.bt
                / physics_variables.plascur,
                "OP ",
            )

            po.ovarrf(
                self.outfile,
                "Normalised total beta",
                " ",
                physics_variables.normalised_total_beta,
                "OP ",
            )

            normalised_toroidal_beta = (
                physics_variables.normalised_total_beta
                * (physics_variables.btot / physics_variables.bt) ** 2
            )
            po.ovarrf(
                self.outfile,
                "Normalised toroidal beta",
                "(normalised_toroidal_beta)",
                normalised_toroidal_beta,
                "OP ",
            )

        if physics_variables.iculbl == 0:
            po.ovarrf(
                self.outfile,
                "Limit on total beta",
                "(betalim)",
                physics_variables.betalim,
                "OP ",
            )
        elif physics_variables.iculbl == 1:
            po.ovarrf(
                self.outfile,
                "Limit on thermal beta",
                "(betalim)",
                physics_variables.betalim,
                "OP ",
            )
        else:
            po.ovarrf(
                self.outfile,
                "Limit on thermal + NB beta",
                "(betalim)",
                physics_variables.betalim,
                "OP ",
            )

        po.ovarre(
            self.outfile,
            "Plasma thermal energy (J)",
            " ",
            1.5e0
            * betath
            * physics_variables.btot
            * physics_variables.btot
            / (2.0e0 * constants.rmu0)
            * physics_variables.vol,
            "OP ",
        )

        po.ovarre(
            self.outfile,
            "Total plasma internal energy (J)",
            "(total_plasma_internal_energy)",
            physics_module.total_plasma_internal_energy,
            "OP ",
        )

        po.osubhd(self.outfile, "Temperature and Density (volume averaged) :")
        po.ovarrf(
            self.outfile, "Electron temperature (keV)", "(te)", physics_variables.te
        )
        po.ovarrf(
            self.outfile,
            "Electron temperature on axis (keV)",
            "(te0)",
            physics_variables.te0,
            "OP ",
        )
        po.ovarrf(self.outfile, "Ion temperature (keV)", "(ti)", physics_variables.ti)
        po.ovarrf(
            self.outfile,
            "Ion temperature on axis (keV)",
            "(ti0)",
            physics_variables.ti0,
            "OP ",
        )
        po.ovarrf(
            self.outfile,
            "Electron temp., density weighted (keV)",
            "(ten)",
            physics_variables.ten,
            "OP ",
        )
        po.ovarre(
            self.outfile, "Electron density (/m3)", "(dene)", physics_variables.dene
        )
        po.ovarre(
            self.outfile,
            "Electron density on axis (/m3)",
            "(ne0)",
            physics_variables.ne0,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Line-averaged electron density (/m3)",
            "(dnla)",
            physics_variables.dnla,
            "OP ",
        )

        if stellarator_variables.istell == 0:
            po.ovarre(
                self.outfile,
                "Line-averaged electron density / Greenwald density",
                "(dnla_gw)",
                physics_variables.dnla / physics_variables.dlimit[6],
                "OP ",
            )

        po.ovarre(
            self.outfile,
            "Ion density (/m3)",
            "(dnitot)",
            physics_variables.dnitot,
            "OP ",
        )
        po.ovarre(
            self.outfile, "Fuel density (/m3)", "(deni)", physics_variables.deni, "OP "
        )
        po.ovarre(
            self.outfile,
            "Total impurity density with Z > 2 (no He) (/m3)",
            "(dnz)",
            physics_variables.dnz,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Helium ion density (thermalised ions only) (/m3)",
            "(dnalp)",
            physics_variables.dnalp,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Proton density (/m3)",
            "(dnprot)",
            physics_variables.dnprot,
            "OP ",
        )
        if physics_variables.protium > 1.0e-10:
            po.ovarre(
                self.outfile,
                "Seeded protium density / electron density",
                "(protium)",
                physics_variables.protium,
            )

        po.ovarre(
            self.outfile,
            "Hot beam density (/m3)",
            "(dnbeam)",
            physics_variables.dnbeam,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Density limit from scaling (/m3)",
            "(dnelimt)",
            physics_variables.dnelimt,
            "OP ",
        )
        if (numerics.ioptimz > 0) and (numerics.active_constraints[4]):
            po.ovarre(
                self.outfile,
                "Density limit (enforced) (/m3)",
                "(boundu(9)*dnelimt)",
                numerics.boundu[8] * physics_variables.dnelimt,
                "OP ",
            )

        po.ovarre(
            self.outfile,
            "Helium ion density (thermalised ions only) / electron density",
            "(ralpne)",
            physics_variables.ralpne,
        )
        po.oblnkl(self.outfile)

        po.ocmmnt(self.outfile, "Impurities")
        po.oblnkl(self.outfile)
        po.ocmmnt(self.outfile, "Plasma ion densities / electron density:")

        for imp in range(impurity_radiation_module.nimp):
            # MDK Update fimp, as this will make the ITV output work correctly.
            impurity_radiation_module.fimp[
                imp
            ] = impurity_radiation_module.impurity_arr_frac[imp]
            str1 = (
                f2py_compatible_to_string(
                    impurity_radiation_module.impurity_arr_label[imp]
                )
                + " concentration"
            )
            str2 = f"(fimp({imp + 1:02}))"
            # MDK Add output flag for H which is calculated.
            if imp == 0:
                po.ovarre(
                    self.outfile, str1, str2, impurity_radiation_module.fimp[imp], "OP "
                )
            else:
                po.ovarre(self.outfile, str1, str2, impurity_radiation_module.fimp[imp])

        po.ovarre(
            self.outfile,
            "Average mass of all ions (amu)",
            "(aion)",
            physics_variables.aion,
            "OP ",
        )
        # MDK Say which impurity is varied, if iteration variable fimpvar (102) is turned on
        # if (any(ixc == 102)) :
        #   call ovarst(self.outfile,'Impurity used as an iteration variable' , '', '"' // impurity_arr(impvar)%label // '"')
        #   po.ovarre(self.outfile,'Fractional density of variable impurity (ion / electron density)','(fimpvar)',fimpvar)
        #
        po.oblnkl(self.outfile)
        po.ovarrf(
            self.outfile, "Effective charge", "(zeff)", physics_variables.zeff, "OP "
        )

        # Issue #487.  No idea what zeffai is.
        # I haven't removed it as it is used in subroutine rether,
        #  (routine to find the equilibration power between the ions and electrons)
        # po.ovarrf(self.outfile,'Mass weighted effective charge','(zeffai)',zeffai, 'OP ')

        po.ovarrf(
            self.outfile, "Density profile factor", "(alphan)", physics_variables.alphan
        )
        po.ovarin(
            self.outfile,
            "Plasma profile model",
            "(ipedestal)",
            physics_variables.ipedestal,
        )

        if physics_variables.ipedestal >= 1:
            if physics_variables.ne0 < physics_variables.neped:
                error_handling.report_error(213)

            po.ocmmnt(self.outfile, "Pedestal profiles are used.")
            po.ovarrf(
                self.outfile,
                "Density pedestal r/a location",
                "(rhopedn)",
                physics_variables.rhopedn,
            )
            if physics_variables.fgwped >= 0e0:
                po.ovarre(
                    self.outfile,
                    "Electron density pedestal height (/m3)",
                    "(neped)",
                    physics_variables.neped,
                    "OP ",
                )
            else:
                po.ovarre(
                    self.outfile,
                    "Electron density pedestal height (/m3)",
                    "(neped)",
                    physics_variables.neped,
                )

            # This code is ODD# Don't change it# No explanation why fgwped and physics_variables.fgwsep
            # must be assigned to their exisiting values#
            fgwped_out = physics_variables.neped / physics_variables.dlimit[6]
            fgwsep_out = physics_variables.nesep / physics_variables.dlimit[6]
            if physics_variables.fgwped >= 0e0:
                physics_variables.fgwped = (
                    physics_variables.neped / physics_variables.dlimit[6]
                )
            if physics_variables.fgwsep >= 0e0:
                physics_variables.fgwsep = (
                    physics_variables.nesep / physics_variables.dlimit[6]
                )

            po.ovarre(
                self.outfile,
                "Electron density at pedestal / nGW",
                "(fgwped_out)",
                fgwped_out,
            )
            po.ovarrf(
                self.outfile,
                "Temperature pedestal r/a location",
                "(rhopedt)",
                physics_variables.rhopedt,
            )

            po.ovarrf(
                self.outfile,
                "Electron temp. pedestal height (keV)",
                "(teped)",
                physics_variables.teped,
            )
            if any(numerics.icc == 78):
                po.ovarrf(
                    self.outfile,
                    "Electron temp. at separatrix (keV)",
                    "(tesep)",
                    physics_variables.tesep,
                    "OP ",
                )
            else:
                po.ovarrf(
                    self.outfile,
                    "Electron temp. at separatrix (keV)",
                    "(tesep)",
                    physics_variables.tesep,
                )

            po.ovarre(
                self.outfile,
                "Electron density at separatrix (/m3)",
                "(nesep)",
                physics_variables.nesep,
            )
            po.ovarre(
                self.outfile,
                "Electron density at separatrix / nGW",
                "(fgwsep_out)",
                fgwsep_out,
            )

        # Issue 558 - addition of constraint 76 to limit the value of nesep, in proportion with the ballooning parameter and Greenwald density
        if any(numerics.icc == 76):
            po.ovarre(
                self.outfile,
                "Critical ballooning parameter value",
                "(alpha_crit)",
                physics_variables.alpha_crit,
            )
            po.ovarre(
                self.outfile,
                "Critical electron density at separatrix (/m3)",
                "(nesep_crit)",
                physics_variables.nesep_crit,
            )

        po.ovarrf(
            self.outfile,
            "Temperature profile index",
            "(alphat)",
            physics_variables.alphat,
        )
        po.ovarrf(
            self.outfile,
            "Temperature profile index beta",
            "(tbeta)",
            physics_variables.tbeta,
        )

        if stellarator_variables.istell == 0:
            po.osubhd(self.outfile, "Density Limit using different models :")
            po.ovarre(
                self.outfile,
                "Old ASDEX model",
                "(dlimit(1))",
                physics_variables.dlimit[0],
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Borrass ITER model I",
                "(dlimit(2))",
                physics_variables.dlimit[1],
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Borrass ITER model II",
                "(dlimit(3))",
                physics_variables.dlimit[2],
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "JET edge radiation model",
                "(dlimit(4))",
                physics_variables.dlimit[3],
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "JET simplified model",
                "(dlimit(5))",
                physics_variables.dlimit[4],
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Hugill-Murakami Mq model",
                "(dlimit(6))",
                physics_variables.dlimit[5],
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Greenwald model",
                "(dlimit(7))",
                physics_variables.dlimit[6],
                "OP ",
            )

        po.osubhd(self.outfile, "Fuel Constituents :")
        po.ovarrf(
            self.outfile, "Deuterium fuel fraction", "(fdeut)", physics_variables.fdeut
        )
        po.ovarrf(
            self.outfile, "Tritium fuel fraction", "(ftrit)", physics_variables.ftrit
        )
        if physics_variables.fhe3 > 1.0e-3:
            po.ovarrf(
                self.outfile, "3-Helium fuel fraction", "(fhe3)", physics_variables.fhe3
            )

        po.osubhd(self.outfile, "Fusion Power :")
        po.ovarre(
            self.outfile,
            "Total fusion power (MW)",
            "(powfmw)",
            physics_variables.powfmw,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            " =    D-T fusion power (MW)",
            "(pdt)",
            physics_variables.pdt,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "  +   D-D fusion power (MW)",
            "(pdd)",
            physics_variables.pdd,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "  + D-He3 fusion power (MW)",
            "(pdhe3)",
            physics_variables.pdhe3,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Alpha power: total (MW)",
            "(palpmw)",
            physics_variables.palpmw,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Alpha power: beam-plasma (MW)",
            "(palpnb)",
            physics_variables.palpnb,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Neutron power (MW)",
            "(pneutmw)",
            physics_variables.pneutmw,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Charged particle power (excluding alphas) (MW)",
            "(pchargemw)",
            physics_variables.pchargemw,
            "OP ",
        )
        tot_power_plasma = (
            physics_variables.falpha * physics_variables.palpmw
            + physics_variables.pchargemw
            + physics_variables.pohmmw
            + current_drive_variables.pinjmw
        )
        po.ovarre(
            self.outfile,
            "Total power deposited in plasma (MW)",
            "(tot_power_plasma)",
            tot_power_plasma,
            "OP ",
        )
        # po.ovarre(self.outfile,'Total power deposited in plasma (MW)','()',falpha*palpmw+pchargemw+pohmmw+pinjmw, 'OP ')

        po.osubhd(self.outfile, "Radiation Power (excluding SOL):")
        po.ovarre(
            self.outfile,
            "Synchrotron radiation power (MW)",
            "(psyncpv*vol)",
            physics_variables.psyncpv * physics_variables.vol,
            "OP ",
        )
        po.ovarrf(
            self.outfile,
            "Synchrotron wall reflectivity factor",
            "(ssync)",
            physics_variables.ssync,
        )
        po.ovarre(
            self.outfile,
            "Normalised minor radius defining 'core'",
            "(coreradius)",
            impurity_radiation_module.coreradius,
        )
        po.ovarre(
            self.outfile,
            "Fraction of core radiation subtracted from P_L",
            "(coreradiationfraction)",
            impurity_radiation_module.coreradiationfraction,
        )
        po.ovarre(
            self.outfile,
            "Radiation power from inner zone (MW)",
            "(pinnerzoneradmw)",
            physics_variables.pinnerzoneradmw,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Radiation power from outer zone (MW)",
            "(pouterzoneradmw)",
            physics_variables.pouterzoneradmw,
            "OP ",
        )

        if stellarator_variables.istell != 0:
            po.ovarre(
                self.outfile,
                "SOL radiation power as imposed by f_rad (MW)",
                "(psolradmw)",
                physics_variables.psolradmw,
                "OP ",
            )

        po.ovarre(
            self.outfile,
            "Total radiation power from inside LCFS (MW)",
            "(pradmw)",
            physics_variables.pradmw,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "LCFS radiation fraction = total radiation in LCFS / total power deposited in plasma",
            "(rad_fraction_LCFS)",
            physics_module.rad_fraction_lcfs,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Nominal mean radiation load on inside surface of reactor (MW/m2)",
            "(photon_wall)",
            physics_variables.photon_wall,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Peaking factor for radiation wall load",
            "(peakfactrad)",
            constraint_variables.peakfactrad,
            "IP ",
        )
        po.ovarre(
            self.outfile,
            "Maximum permitted radiation wall load (MW/m^2)",
            "(maxradwallload)",
            constraint_variables.maxradwallload,
            "IP ",
        )
        po.ovarre(
            self.outfile,
            "Peak radiation wall load (MW/m^2)",
            "(peakradwallload)",
            constraint_variables.peakradwallload,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Fast alpha particle power incident on the first wall (MW)",
            "(palpfwmw)",
            physics_variables.palpfwmw,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Nominal mean neutron load on inside surface of reactor (MW/m2)",
            "(wallmw)",
            physics_variables.wallmw,
            "OP ",
        )

        if stellarator_variables.istell == 0:
            po.oblnkl(self.outfile)
            po.ovarre(
                self.outfile,
                "Power incident on the divertor targets (MW)",
                "(ptarmw)",
                physics_module.ptarmw,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Fraction of power to the lower divertor",
                "(ftar)",
                physics_variables.ftar,
                "IP ",
            )
            po.ovarre(
                self.outfile,
                "Outboard side heat flux decay length (m)",
                "(lambdaio)",
                physics_module.lambdaio,
                "OP ",
            )
            if physics_variables.idivrt == 2:
                po.ovarre(
                    self.outfile,
                    "Midplane seperation of the two magnetic closed flux surfaces (m)",
                    "(drsep)",
                    physics_module.drsep,
                    "OP ",
                )

            po.ovarre(
                self.outfile,
                "Fraction of power on the inner targets",
                "(fio)",
                physics_module.fio,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Fraction of power incident on the lower inner target",
                "(fLI)",
                physics_module.fli,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Fraction of power incident on the lower outer target",
                "(fLO)",
                physics_module.flo,
                "OP ",
            )
            if physics_variables.idivrt == 2:
                po.ovarre(
                    self.outfile,
                    "Fraction of power incident on the upper inner target",
                    "(fUI)",
                    physics_module.fui,
                    "OP ",
                )
                po.ovarre(
                    self.outfile,
                    "Fraction of power incident on the upper outer target",
                    "(fUO)",
                    physics_module.fuo,
                    "OP ",
                )

            po.ovarre(
                self.outfile,
                "Power incident on the lower inner target (MW)",
                "(pLImw)",
                physics_module.plimw,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Power incident on the lower outer target (MW)",
                "(pLOmw)",
                physics_module.plomw,
                "OP ",
            )
            if physics_variables.idivrt == 2:
                po.ovarre(
                    self.outfile,
                    "Power incident on the upper innner target (MW)",
                    "(pUImw)",
                    physics_module.puimw,
                    "OP ",
                )
                po.ovarre(
                    self.outfile,
                    "Power incident on the upper outer target (MW)",
                    "(pUOmw)",
                    physics_module.puomw,
                    "OP ",
                )

        po.oblnkl(self.outfile)
        po.ovarre(
            self.outfile,
            "Ohmic heating power (MW)",
            "(pohmmw)",
            physics_variables.pohmmw,
            "OP ",
        )
        po.ovarrf(
            self.outfile,
            "Fraction of alpha power deposited in plasma",
            "(falpha)",
            physics_variables.falpha,
            "IP",
        )
        po.ovarrf(
            self.outfile,
            "Fraction of alpha power to electrons",
            "(falpe)",
            physics_variables.falpe,
            "OP ",
        )
        po.ovarrf(
            self.outfile,
            "Fraction of alpha power to ions",
            "(falpi)",
            physics_variables.falpi,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Ion transport (MW)",
            "(ptrimw)",
            physics_variables.ptrimw,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Electron transport (MW)",
            "(ptremw)",
            physics_variables.ptremw,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Injection power to ions (MW)",
            "(pinjimw)",
            current_drive_variables.pinjimw,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Injection power to electrons (MW)",
            "(pinjemw)",
            current_drive_variables.pinjemw,
            "OP ",
        )
        if physics_variables.ignite == 1:
            po.ocmmnt(self.outfile, "  (Injected power only used for start-up phase)")

        po.ovarin(
            self.outfile,
            "Ignited plasma switch (0=not ignited, 1=ignited)",
            "(ignite)",
            physics_variables.ignite,
        )

        po.oblnkl(self.outfile)
        po.ovarre(
            self.outfile,
            "Power into divertor zone via charged particles (MW)",
            "(pdivt)",
            physics_variables.pdivt,
            "OP ",
        )

        if physics_variables.pdivt <= 0.001e0:
            error_handling.fdiags[0] = physics_variables.pdivt
            error_handling.report_error(87)
            po.oblnkl(self.outfile)
            po.ocmmnt(
                self.outfile, "  BEWARE: possible problem with high radiation power"
            )
            po.ocmmnt(
                self.outfile, "          Power into divertor zone is unrealistic;"
            )
            po.ocmmnt(self.outfile, "          divertor calculations will be nonsense#")
            po.ocmmnt(
                self.outfile, "  Set constraint 17 (Radiation fraction upper limit)."
            )
            po.oblnkl(self.outfile)

        if physics_variables.idivrt == 2:
            # Double null divertor configuration
            po.ovarre(
                self.outfile,
                "Pdivt / R ratio (MW/m) (On peak divertor)",
                "(pdivmax/physics_variables.rmajor)",
                physics_variables.pdivmax / physics_variables.rmajor,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Pdivt Bt / qAR ratio (MWT/m) (On peak divertor)",
                "(pdivmaxbt/qar)",
                (
                    (physics_variables.pdivmax * physics_variables.bt)
                    / (
                        physics_variables.q95
                        * physics_variables.aspect
                        * physics_variables.rmajor
                    )
                ),
                "OP ",
            )
        else:
            # Single null divertor configuration
            po.ovarre(
                self.outfile,
                "Psep / R ratio (MW/m)",
                "(pdivt/rmajor)",
                physics_variables.pdivt / physics_variables.rmajor,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Psep Bt / qAR ratio (MWT/m)",
                "(pdivtbt/qar)",
                (
                    (physics_variables.pdivt * physics_variables.bt)
                    / (
                        physics_variables.q95
                        * physics_variables.aspect
                        * physics_variables.rmajor
                    )
                ),
                "OP ",
            )

        if stellarator_variables.istell == 0:
            po.osubhd(self.outfile, "H-mode Power Threshold Scalings :")

            po.ovarre(
                self.outfile,
                "ITER 1996 scaling: nominal (MW)",
                "(pthrmw(1))",
                physics_variables.pthrmw[0],
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "ITER 1996 scaling: upper bound (MW)",
                "(pthrmw(2))",
                physics_variables.pthrmw[1],
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "ITER 1996 scaling: lower bound (MW)",
                "(pthrmw(3))",
                physics_variables.pthrmw[2],
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "ITER 1997 scaling (1) (MW)",
                "(pthrmw(4))",
                physics_variables.pthrmw[3],
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "ITER 1997 scaling (2) (MW)",
                "(pthrmw(5))",
                physics_variables.pthrmw[4],
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Martin 2008 scaling: nominal (MW)",
                "(pthrmw(6))",
                physics_variables.pthrmw[5],
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Martin 2008 scaling: 95% upper bound (MW)",
                "(pthrmw(7))",
                physics_variables.pthrmw[6],
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Martin 2008 scaling: 95% lower bound (MW)",
                "(pthrmw(8))",
                physics_variables.pthrmw[7],
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Snipes 2000 scaling: nominal (MW)",
                "(pthrmw(9))",
                physics_variables.pthrmw[8],
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Snipes 2000 scaling: upper bound (MW)",
                "(pthrmw(10))",
                physics_variables.pthrmw[9],
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Snipes 2000 scaling: lower bound (MW)",
                "(pthrmw(11))",
                physics_variables.pthrmw[10],
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Snipes 2000 scaling (closed divertor): nominal (MW)",
                "(pthrmw(12))",
                physics_variables.pthrmw[11],
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Snipes 2000 scaling (closed divertor): upper bound (MW)",
                "(pthrmw(13))",
                physics_variables.pthrmw[12],
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Snipes 2000 scaling (closed divertor): lower bound (MW)",
                "(pthrmw(14))",
                physics_variables.pthrmw[13],
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Hubbard 2012 L-I threshold - nominal (MW)",
                "(pthrmw(15))",
                physics_variables.pthrmw[14],
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Hubbard 2012 L-I threshold - lower bound (MW)",
                "(pthrmw(16))",
                physics_variables.pthrmw[15],
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Hubbard 2012 L-I threshold - upper bound (MW)",
                "(pthrmw(17))",
                physics_variables.pthrmw[16],
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Hubbard 2017 L-I threshold",
                "(pthrmw(18))",
                physics_variables.pthrmw[17],
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Martin 2008 aspect ratio corrected scaling: nominal (MW)",
                "(pthrmw(19))",
                physics_variables.pthrmw[18],
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Martin 2008 aspect ratio corrected scaling: 95% upper bound (MW)",
                "(pthrmw(20))",
                physics_variables.pthrmw[19],
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Martin 2008 aspect ratio corrected scaling: 95% lower bound (MW)",
                "(pthrmw(21))",
                physics_variables.pthrmw[20],
                "OP ",
            )
            po.oblnkl(self.outfile)
            if physics_variables.ilhthresh in [9, 10, 11]:
                if (physics_variables.bt < 0.78e0) or (physics_variables.bt > 7.94e0):
                    po.ocmmnt(
                        self.outfile,
                        "(physics_variables.bt outside Snipes 2000 fitted range)",
                    )
                    error_handling.report_error(201)

                if (physics_variables.rminor < 0.15e0) or (
                    physics_variables.rminor > 1.15e0
                ):
                    po.ocmmnt(self.outfile, "(rminor outside Snipes 2000 fitted range)")
                    error_handling.report_error(202)

                if (physics_variables.rmajor < 0.55e0) or (
                    physics_variables.rmajor > 3.37e0
                ):
                    po.ocmmnt(
                        self.outfile,
                        "(physics_variables.rmajor outside Snipes 2000 fitted range)",
                    )
                    error_handling.report_error(203)

                if (physics_variables.dnla < 0.09e20) or (
                    physics_variables.dnla > 3.16e20
                ):
                    po.ocmmnt(
                        self.outfile,
                        "(physics_variables.dnla outside Snipes 2000 fitted range)",
                    )
                    error_handling.report_error(204)

                if (physics_variables.kappa < 1.0e0) or (
                    physics_variables.kappa > 2.04e0
                ):
                    po.ocmmnt(
                        self.outfile,
                        "(physics_variables.kappa outside Snipes 2000 fitted range)",
                    )
                    error_handling.report_error(205)

                if (physics_variables.triang < 0.07e0) or (
                    physics_variables.triang > 0.74e0
                ):
                    po.ocmmnt(self.outfile, "(triang outside Snipes 2000 fitted range)")
                    error_handling.report_error(206)

            po.oblnkl(self.outfile)

            if physics_variables.ilhthresh in [12, 13, 14]:
                po.ocmmnt(
                    self.outfile,
                    "(L-H threshold for closed divertor only. Limited data used in Snipes fit)",
                )
                po.oblnkl(self.outfile)
                error_handling.report_error(207)

            if (numerics.ioptimz > 0) and (numerics.active_constraints[14]):
                po.ovarre(
                    self.outfile,
                    "L-H threshold power (enforced) (MW)",
                    "(boundl(103)*plhthresh)",
                    numerics.boundl[102] * physics_variables.plhthresh,
                    "OP ",
                )
                po.ovarre(
                    self.outfile,
                    "L-H threshold power (MW)",
                    "(plhthresh)",
                    physics_variables.plhthresh,
                    "OP ",
                )
            else:
                po.ovarre(
                    self.outfile,
                    "L-H threshold power (NOT enforced) (MW)",
                    "(plhthresh)",
                    physics_variables.plhthresh,
                    "OP ",
                )

        po.osubhd(self.outfile, "Confinement :")

        if physics_variables.ignite == 1:
            po.ocmmnt(
                self.outfile,
                "Device is assumed to be ignited for the calculation of confinement time",
            )
            po.oblnkl(self.outfile)

        tauelaw = f2py_compatible_to_string(
            physics_variables.tauscl[physics_variables.isc - 1]
        )

        po.ocmmnt(
            self.outfile,
            f"Confinement scaling law: {tauelaw}",
        )

        po.ovarst(
            self.outfile,
            "Confinement scaling law",
            "(tauelaw)",
            f'"{tauelaw.strip().split(" ")[0]}"',
        )

        po.ovarrf(
            self.outfile, "Confinement H factor", "(hfact)", physics_variables.hfact
        )
        po.ovarrf(
            self.outfile,
            "Global thermal energy confinement time (s)",
            "(taueff)",
            physics_variables.taueff,
            "OP ",
        )
        po.ovarrf(
            self.outfile,
            "Ion energy confinement time (s)",
            "(tauei)",
            physics_variables.tauei,
            "OP ",
        )
        po.ovarrf(
            self.outfile,
            "Electron energy confinement time (s)",
            "(tauee)",
            physics_variables.tauee,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "n.tau = Volume-average electron density x Energy confinement time (s/m3)",
            "(dntau)",
            physics_variables.dntau,
            "OP ",
        )
        po.ocmmnt(
            self.outfile,
            "Triple product = Vol-av electron density x Vol-av electron temp x Energy confinement time:",
        )
        po.ovarre(
            self.outfile,
            "Triple product  (keV s/m3)",
            "(dntau*te)",
            physics_variables.dntau * physics_variables.te,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Transport loss power assumed in scaling law (MW)",
            "(powerht)",
            physics_variables.powerht,
            "OP ",
        )
        po.ovarin(
            self.outfile,
            "Switch for radiation loss term usage in power balance",
            "(iradloss)",
            physics_variables.iradloss,
        )
        if physics_variables.iradloss == 0:
            po.ovarre(
                self.outfile,
                "Radiation power subtracted from plasma power balance (MW)",
                "",
                physics_variables.pradmw,
                "OP ",
            )
            po.ocmmnt(self.outfile, "  (Radiation correction is total radiation power)")
        elif physics_variables.iradloss == 1:
            po.ovarre(
                self.outfile,
                "Radiation power subtracted from plasma power balance (MW)",
                "",
                physics_variables.pinnerzoneradmw,
                "OP ",
            )
            po.ocmmnt(self.outfile, "  (Radiation correction is core radiation power)")
        else:
            po.ovarre(
                self.outfile,
                "Radiation power subtracted from plasma power balance (MW)",
                "",
                0.0e0,
            )
            po.ocmmnt(self.outfile, "  (No radiation correction applied)")
        if physics_variables.iradloss == 1:
            po.ovarrf(
                self.outfile,
                "H* non-radiation corrected",
                "(hstar)",
                physics_variables.hfact
                * (
                    physics_variables.powerht
                    / (
                        physics_variables.powerht
                        + physics_variables.psyncpv
                        + physics_variables.pinnerzoneradmw
                    )
                )
                ** 0.31,
                "OP",
            )
        elif physics_variables.iradloss == 0:
            po.ovarrf(
                self.outfile,
                "H* non-radiation corrected",
                "(hstar)",
                physics_variables.hfact
                * (
                    physics_variables.powerht
                    / (
                        physics_variables.powerht
                        + physics_variables.pradpv * physics_variables.vol
                    )
                )
                ** 0.31,
                "OP",
            )
        elif physics_variables.iradloss == 2:
            po.ovarrf(
                self.outfile,
                "H* non-radiation corrected",
                "(hstar)",
                physics_variables.hfact,
                "OP",
            )
        po.ocmmnt(self.outfile, "  (H* assumes IPB98(y,2), ELMy H-mode scaling)")
        po.ovarrf(
            self.outfile,
            "Alpha particle confinement time (s)",
            "(taup)",
            physics_variables.taup,
            "OP ",
        )
        # Note alpha confinement time is no longer equal to fuel particle confinement time.
        po.ovarrf(
            self.outfile,
            "Alpha particle/energy confinement time ratio",
            "(taup/taueff)",
            physics_variables.taup / physics_variables.taueff,
            "OP ",
        )
        po.ovarrf(
            self.outfile,
            "Lower limit on taup/taueff",
            "(taulimit)",
            constraint_variables.taulimit,
        )
        po.ovarrf(
            self.outfile,
            "Total energy confinement time including radiation loss (s)",
            "(total_energy_conf_time)",
            physics_module.total_energy_conf_time,
            "OP ",
        )
        po.ocmmnt(
            self.outfile,
            "  (= stored energy including fast particles / loss power including radiation",
        )

        if stellarator_variables.istell == 0:
            # Issues 363 Output dimensionless plasma parameters MDK
            po.osubhd(self.outfile, "Dimensionless plasma parameters")
            po.ocmmnt(self.outfile, "For definitions see")
            po.ocmmnt(
                self.outfile,
                "Recent progress on the development and analysis of the ITPA global H-mode confinement database",
            )
            po.ocmmnt(
                self.outfile,
                "D.C. McDonald et al, 2007 Nuclear Fusion v47, 147. (nu_star missing 1/mu0)",
            )
            po.ovarre(
                self.outfile,
                "Normalized plasma pressure beta as defined by McDonald et al",
                "(beta_mcdonald)",
                physics_module.beta_mcdonald,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Normalized ion Larmor radius",
                "(rho_star)",
                physics_module.rho_star,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Normalized collisionality",
                "(nu_star)",
                physics_module.nu_star,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Volume measure of elongation",
                "(kappaa_IPB)",
                physics_variables.kappaa_ipb,
                "OP ",
            )

            po.osubhd(self.outfile, "Plasma Volt-second Requirements :")
            po.ovarre(
                self.outfile,
                "Total volt-second requirement (Wb)",
                "(vsstt)",
                physics_variables.vsstt,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Inductive volt-seconds (Wb)",
                "(vsind)",
                physics_variables.vsind,
                "OP ",
            )
            po.ovarrf(
                self.outfile, "Ejima coefficient", "(gamma)", physics_variables.gamma
            )
            po.ovarre(
                self.outfile,
                "Start-up resistive (Wb)",
                "(vsres)",
                physics_variables.vsres,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Flat-top resistive (Wb)",
                "(vsbrn)",
                physics_variables.vsbrn,
                "OP ",
            )

            po.ovarrf(
                self.outfile,
                "Bootstrap current fraction multiplier",
                "(cboot)",
                current_drive_variables.cboot,
            )
            po.ovarrf(
                self.outfile,
                "Bootstrap fraction (ITER 1989)",
                "(bscf_iter89)",
                current_drive_variables.bscf_iter89,
                "OP ",
            )

            po.ovarrf(
                self.outfile,
                "Bootstrap fraction (Sauter et al)",
                "(bscf_sauter)",
                current_drive_variables.bscf_sauter,
                "OP ",
            )

            po.ovarrf(
                self.outfile,
                "Bootstrap fraction (Nevins et al)",
                "(bscf_nevins)",
                current_drive_variables.bscf_nevins,
                "OP ",
            )
            po.ovarrf(
                self.outfile,
                "Bootstrap fraction (Wilson)",
                "(bscf_wilson)",
                current_drive_variables.bscf_wilson,
                "OP ",
            )

            po.ovarrf(
                self.outfile,
                "Bootstrap fraction (Sakai)",
                "(bscf_sakai)",
                current_drive_variables.bscf_sakai,
                "OP ",
            )
            po.ovarrf(
                self.outfile,
                "Diamagnetic fraction (Hender)",
                "(diacf_hender)",
                current_drive_variables.diacf_hender,
                "OP ",
            )
            po.ovarrf(
                self.outfile,
                "Diamagnetic fraction (SCENE)",
                "(diacf_scene)",
                current_drive_variables.diacf_scene,
                "OP ",
            )
            po.ovarrf(
                self.outfile,
                "Pfirsch-Schlueter fraction (SCENE)",
                "(pscf_scene)",
                current_drive_variables.pscf_scene,
                "OP ",
            )
            # Error to catch if bootstap fraction limit has been enforced
            if physics_module.err242 == 1:
                error_handling.report_error(242)

            # Error to catch if self-driven current fraction limit has been enforced
            if physics_module.err243 == 1:
                error_handling.report_error(243)

            if current_drive_variables.bscfmax < 0.0e0:
                po.ocmmnt(
                    self.outfile, "  (User-specified bootstrap current fraction used)"
                )
            elif physics_variables.ibss == 1:
                po.ocmmnt(
                    self.outfile, "  (ITER 1989 bootstrap current fraction model used)"
                )
            elif physics_variables.ibss == 2:
                po.ocmmnt(
                    self.outfile,
                    "  (Nevins et al bootstrap current fraction model used)",
                )
            elif physics_variables.ibss == 3:
                po.ocmmnt(
                    self.outfile, "  (Wilson bootstrap current fraction model used)"
                )
            elif physics_variables.ibss == 4:
                po.ocmmnt(
                    self.outfile,
                    "  (Sauter et al bootstrap current fraction model used)",
                )
            elif physics_variables.ibss == 5:
                po.ocmmnt(
                    self.outfile,
                    "  (Sakai et al bootstrap current fraction model used)",
                )

            if physics_variables.idia == 0:
                po.ocmmnt(
                    self.outfile, "  (Diamagnetic current fraction not calculated)"
                )
                # Error to show if diamagnetic current is above 1% but not used
                if current_drive_variables.diacf_scene > 0.01e0:
                    error_handling.report_error(244)

            elif physics_variables.idia == 1:
                po.ocmmnt(
                    self.outfile, "  (Hender diamagnetic current fraction scaling used)"
                )
            elif physics_variables.idia == 2:
                po.ocmmnt(
                    self.outfile, "  (SCENE diamagnetic current fraction scaling used)"
                )

            if physics_variables.ips == 0:
                po.ocmmnt(
                    self.outfile, "  Pfirsch-Schluter current fraction not calculated"
                )
            elif physics_variables.ips == 1:
                po.ocmmnt(
                    self.outfile,
                    "  (SCENE Pfirsch-Schluter current fraction scaling used)",
                )

            po.ovarrf(
                self.outfile,
                "Bootstrap fraction (enforced)",
                "(bootipf.)",
                current_drive_variables.bootipf,
                "OP ",
            )
            po.ovarrf(
                self.outfile,
                "Diamagnetic fraction (enforced)",
                "(diaipf.)",
                current_drive_variables.diaipf,
                "OP ",
            )
            po.ovarrf(
                self.outfile,
                "Pfirsch-Schlueter fraction (enforced)",
                "(psipf.)",
                current_drive_variables.psipf,
                "OP ",
            )

            po.ovarre(
                self.outfile,
                "Loop voltage during burn (V)",
                "(vburn)",
                physics_variables.plascur
                * physics_variables.rplas
                * physics_variables.facoh,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Plasma resistance (ohm)",
                "(rplas)",
                physics_variables.rplas,
                "OP ",
            )

            po.ovarre(
                self.outfile,
                "Resistive diffusion time (s)",
                "(res_time)",
                physics_variables.res_time,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Plasma inductance (H)",
                "(rlp)",
                physics_variables.rlp,
                "OP ",
            )
            po.ovarrf(
                self.outfile,
                "Coefficient for sawtooth effects on burn V-s requirement",
                "(csawth)",
                physics_variables.csawth,
            )

        po.osubhd(self.outfile, "Fuelling :")
        po.ovarre(
            self.outfile,
            "Ratio of He and pellet particle confinement times",
            "(tauratio)",
            physics_variables.tauratio,
        )
        po.ovarre(
            self.outfile,
            "Fuelling rate (nucleus-pairs/s)",
            "(qfuel)",
            physics_variables.qfuel,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Fuel burn-up rate (reactions/s)",
            "(rndfuel)",
            physics_variables.rndfuel,
            "OP ",
        )
        po.ovarrf(
            self.outfile,
            "Burn-up fraction",
            "(burnup)",
            physics_variables.burnup,
            "OP ",
        )

        if any(numerics.icc == 78):
            po.osubhd(self.outfile, "Reinke Criterion :")
            po.ovarin(
                self.outfile,
                "index of impurity to be iterated for divertor detachment",
                "(impvardiv)",
                reinke_variables.impvardiv,
            )
            po.ovarre(
                self.outfile,
                "Minimum Impurity fraction from Reinke",
                "(fzmin)",
                reinke_variables.fzmin,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Actual Impurity fraction",
                "(fzactual)",
                reinke_variables.fzactual,
            )

    def igmarcal(self):
        """Routine to calculate ignition margin
        author: P J Knight, CCFE, Culham Science Centre
        outfile   : input integer : Fortran output unit identifier
        This routine calculates the ignition margin at the final point
        with different scalings.
        """

        po.oheadr(self.outfile, "Energy confinement times, and required H-factors :")
        po.ocmmnt(
            self.outfile,
            f"{'':>5}{'scaling law':<25}{'confinement time (s)':<25}H-factor for",
        )
        po.ocmmnt(
            self.outfile,
            f"{'':>34}{'for H = 1':<20}power balance",
        )

        # for iisc in range(32, 48):
        # Put the ITPA value first
        for iisc in [50, 34, 37, 38, 39, 46, 47, 48]:
            (
                physics_variables.kappaa,
                ptrez,
                ptriz,
                taueez,
                taueiz,
                taueffz,
                powerhtz,
            ) = self.pcond(
                physics_variables.afuel,
                physics_variables.palpmw,
                physics_variables.aspect,
                physics_variables.bt,
                physics_variables.dnitot,
                physics_variables.dene,
                physics_variables.dnla,
                physics_variables.eps,
                1.0,
                physics_variables.iinvqd,
                iisc,
                physics_variables.ignite,
                physics_variables.kappa,
                physics_variables.kappa95,
                physics_variables.pchargemw,
                current_drive_variables.pinjmw,
                physics_variables.plascur,
                physics_variables.pcoreradpv,
                physics_variables.rmajor,
                physics_variables.rminor,
                physics_variables.te,
                physics_variables.ten,
                physics_variables.tin,
                physics_variables.q,
                physics_variables.qstar,
                physics_variables.vol,
                physics_variables.xarea,
                physics_variables.zeff,
            )

            physics_variables.hfac[iisc - 1] = self.fhfac(iisc)

            po.ocmmnt(
                self.outfile,
                f"{'':>2}{f2py_compatible_to_string(physics_variables.tauscl[iisc - 1]):<32}"
                f"{taueez:<26.3f}{physics_variables.hfac[iisc - 1]:.3f}",
            )

    @staticmethod
    def bootstrap_fraction_iter89(
        aspect, beta, bt, cboot, plascur, q95, q0, rmajor, vol
    ):
        """Original ITER calculation of bootstrap-driven fraction
        of the plasma current.
        author: P J Knight, CCFE, Culham Science Centre
        aspect  : input real : plasma aspect ratio
        beta    : input real : plasma total beta
        bt      : input real : toroidal field on axis (T)
        cboot   : input real : bootstrap current fraction multiplier
        plascur : input real : plasma current (A)
        q95     : input real : safety factor at 95% surface
        q0      : input real : central safety factor
        rmajor  : input real : plasma major radius (m)
        vol     : input real : plasma volume (m3)
        This routine performs the original ITER calculation of the
        plasma current bootstrap fraction.
        ITER Physics Design Guidelines: 1989 [IPDG89], N. A. Uckan et al,
        ITER Documentation Series No.10, IAEA/ITER/DS/10, IAEA, Vienna, 1990
        """
        xbs = min(10, q95 / q0)
        cbs = cboot * (1.32 - 0.235 * xbs + 0.0185 * xbs**2)
        bpbs = (
            constants.rmu0
            * plascur
            / (2 * np.pi * np.sqrt(vol / (2 * np.pi**2 * rmajor)))
        )
        betapbs = beta * bt**2 / bpbs**2

        if betapbs <= 0.0:  # only possible if beta <= 0.0
            return 0.0
        return cbs * (betapbs / np.sqrt(aspect)) ** 1.3

    @staticmethod
    def bootstrap_fraction_wilson(
        alphaj, alphap, alphat, betpth, q0, qpsi, rmajor, rminor
    ):
        """Bootstrap current fraction from Wilson et al scaling
        author: P J Knight, CCFE, Culham Science Centre
        alphaj  : input real :  current profile index
        alphap  : input real :  pressure profile index
        alphat  : input real :  temperature profile index
        beta    : input real :  total beta
        betpth  : input real :  thermal component of poloidal beta
        q0      : input real :  safety factor on axis
        qpsi    : input real :  edge safety factor
        rmajor  : input real :  major radius (m)
        rminor  : input real :  minor radius (m)
        This function calculates the bootstrap current fraction
        using the numerically fitted algorithm written by Howard Wilson.
        AEA FUS 172: Physics Assessment for the European Reactor Study
        H. R. Wilson, Nuclear Fusion <B>32</B> (1992) 257
        """
        term1 = np.log(0.5)
        term2 = np.log(q0 / qpsi)

        termp = 1.0 - 0.5 ** (1.0 / alphap)
        termt = 1.0 - 0.5 ** (1.0 / alphat)
        termj = 1.0 - 0.5 ** (1.0 / alphaj)

        alfpnw = term1 / np.log(np.log((q0 + (qpsi - q0) * termp) / qpsi) / term2)
        alftnw = term1 / np.log(np.log((q0 + (qpsi - q0) * termt) / qpsi) / term2)
        aj = term1 / np.log(np.log((q0 + (qpsi - q0) * termj) / qpsi) / term2)

        # Crude check for NaN errors or other illegal values...

        if np.isnan(aj) or np.isnan(alfpnw) or np.isnan(alftnw) or aj < 0:
            error_handling.fdiags[0] = aj
            error_handling.fdiags[1] = alfpnw
            error_handling.fdiags[2] = alftnw
            error_handling.fdiags[3] = aj

            error_handling.report_error(76)

        # Ratio of ionic charge to electron charge

        z = 1.0

        # Inverse aspect ratio: r2 = maximum plasma radius, r1 = minimum

        r2 = rmajor + rminor
        r1 = rmajor - rminor
        eps1 = (r2 - r1) / (r2 + r1)

        # Coefficients fitted using least squares techniques
        saj = np.sqrt(aj)

        a = np.array(
            [
                1.41 * (1.0 - 0.28 * saj) * (1.0 + 0.12 / z),
                0.36 * (1.0 - 0.59 * saj) * (1.0 + 0.8 / z),
                -0.27 * (1.0 - 0.47 * saj) * (1.0 + 3.0 / z),
                0.0053 * (1.0 + 5.0 / z),
                -0.93 * (1.0 - 0.34 * saj) * (1.0 + 0.15 / z),
                -0.26 * (1.0 - 0.57 * saj) * (1.0 - 0.27 * z),
                0.064 * (1.0 - 0.6 * aj + 0.15 * aj * aj) * (1.0 + 7.6 / z),
                -0.0011 * (1.0 + 9.0 / z),
                -0.33 * (1.0 - aj + 0.33 * aj * aj),
                -0.26 * (1.0 - 0.87 / saj - 0.16 * aj),
                -0.14 * (1.0 - 1.14 / saj - 0.45 * saj),
                -0.0069,
            ]
        )

        seps1 = np.sqrt(eps1)

        b = np.array(
            [
                1.0,
                alfpnw,
                alftnw,
                alfpnw * alftnw,
                seps1,
                alfpnw * seps1,
                alftnw * seps1,
                alfpnw * alftnw * seps1,
                eps1,
                alfpnw * eps1,
                alftnw * eps1,
                alfpnw * alftnw * eps1,
            ]
        )

        # Empirical bootstrap current fraction
        return seps1 * betpth * (a * b).sum()

    @staticmethod
    def bootstrap_fraction_nevins(
        alphan,
        alphat,
        betat,
        bt,
        dene,
        plascur,
        q95,
        q0,
        rmajor,
        rminor,
        ten,
        zeff,
    ):
        """Bootstrap current fraction from Nevins et al scaling
        author: P J Knight, CCFE, Culham Science Centre
        alphan : input real :  density profile index
        alphat : input real :  temperature profile index
        betat  : input real :  total plasma beta (with respect to the toroidal
        field)
        bt     : input real :  toroidal field on axis (T)
        dene   : input real :  electron density (/m3)
        plascur: input real :  plasma current (A)
        q0     : input real :  central safety factor
        q95    : input real :  safety factor at 95% surface
        rmajor : input real :  plasma major radius (m)
        rminor : input real :  plasma minor radius (m)
        ten    : input real :  density weighted average plasma temperature (keV)
        zeff   : input real :  plasma effective charge
        This function calculates the bootstrap current fraction,
        using the Nevins et al method, 4/11/90.
        """
        # Calculate peak electron beta

        betae0 = (
            physics_variables.ne0
            * physics_variables.te0
            * 1.0e3
            * constants.echarge
            / (bt**2 / (2.0 * constants.rmu0))
        )

        # Call integration routine

        ainteg, _ = integrate.quad(
            lambda y: bsinteg(
                y,
                dene,
                ten,
                bt,
                rminor,
                rmajor,
                zeff,
                alphat,
                alphan,
                q0,
                q95,
                betat,
                constants.echarge,
                constants.rmu0,
            ),
            0,
            0.999,
            epsabs=0.001,
            epsrel=0.001,
        )

        # Calculate bootstrap current and fraction

        aibs = 2.5 * betae0 * rmajor * bt * q95 * ainteg
        return 1.0e6 * aibs / plascur

    @staticmethod
    def bootstrap_fraction_sauter():
        """Bootstrap current fraction from Sauter et al scaling
        author: P J Knight, CCFE, Culham Science Centre
        None
        This function calculates the bootstrap current fraction
        using the Sauter, Angioni and Lin-Liu scaling.
        <P>The code was supplied by Emiliano Fable, IPP Garching
        (private communication).
        O. Sauter, C. Angioni and Y. R. Lin-Liu,
        Physics of Plasmas <B>6</B> (1999) 2834
        O. Sauter, C. Angioni and Y. R. Lin-Liu, (ERRATA)
        Physics of Plasmas <B>9</B> (2002) 5140
        """
        NR = 200

        roa = np.arange(1, NR + 1, step=1) / NR

        rho = np.sqrt(physics_variables.xarea / np.pi) * roa
        sqeps = np.sqrt(roa * (physics_variables.rminor / physics_variables.rmajor))

        ne = 1e-19 * np.vectorize(
            lambda r: profiles_module.nprofile(
                r,
                physics_variables.rhopedn,
                physics_variables.ne0,
                physics_variables.neped,
                physics_variables.nesep,
                physics_variables.alphan,
            )
        )(roa)
        ni = (physics_variables.dnitot / physics_variables.dene) * ne
        tempe = np.vectorize(
            lambda r: profiles_module.tprofile(
                r,
                physics_variables.rhopedt,
                physics_variables.te0,
                physics_variables.teped,
                physics_variables.tesep,
                physics_variables.alphat,
                physics_variables.tbeta,
            )
        )(roa)
        tempi = (physics_variables.ti / physics_variables.te) * tempe

        zef = np.full_like(tempi, physics_variables.zeff)  # Flat Zeff profile assumed

        # mu = 1/safety factor
        # Parabolic q profile assumed
        mu = 1 / (
            physics_variables.q0
            + (physics_variables.q - physics_variables.q0) * roa**2
        )
        amain = np.full_like(mu, physics_variables.afuel)
        zmain = np.full_like(mu, 1.0 + physics_variables.fhe3)

        if ne[NR - 1] == 0.0:
            ne[NR - 1] = 1e-4 * ne[NR - 2]
            ni[NR - 1] = 1e-4 * ni[NR - 2]

        if tempe[NR - 1] == 0.0:
            tempe[NR - 1] = 1e-4 * tempe[NR - 2]
            tempi[NR - 1] = 1e-4 * tempi[NR - 2]

        # Calculate total bootstrap current (MA) by summing along profiles
        # Looping from 2 because dcsa etc should return 0 @ j == 1
        nr_rng = np.arange(2, NR)
        nr_rng_1 = nr_rng - 1

        drho = rho[nr_rng] - rho[nr_rng_1]
        da = 2 * np.pi * rho[nr_rng_1] * drho  # area of annulus
        dlogte_drho = (np.log(tempe[nr_rng]) - np.log(tempe[nr_rng_1])) / drho
        dlogti_drho = (np.log(tempi[nr_rng]) - np.log(tempi[nr_rng_1])) / drho
        dlogne_drho = (np.log(ne[nr_rng]) - np.log(ne[nr_rng_1])) / drho

        jboot = (
            0.5
            * (
                dcsa(
                    nr_rng,
                    NR,
                    physics_variables.rmajor,
                    physics_variables.bt,
                    physics_variables.triang,
                    ne,
                    ni,
                    tempe,
                    tempi,
                    mu,
                    rho,
                    zef,
                    sqeps,
                )
                * dlogne_drho
                + hcsa(
                    nr_rng,
                    NR,
                    physics_variables.rmajor,
                    physics_variables.bt,
                    physics_variables.triang,
                    ne,
                    ni,
                    tempe,
                    tempi,
                    mu,
                    rho,
                    zef,
                    sqeps,
                )
                * dlogte_drho
                + xcsa(
                    nr_rng,
                    NR,
                    physics_variables.rmajor,
                    physics_variables.bt,
                    physics_variables.triang,
                    mu,
                    sqeps,
                    tempi,
                    tempe,
                    amain,
                    zmain,
                    ni,
                    ne,
                    rho,
                    zef,
                )
                * dlogti_drho
            )
            * 1.0e6
            * (
                -physics_variables.bt
                / (0.2 * np.pi * physics_variables.rmajor)
                * rho[nr_rng_1]
                * mu[nr_rng_1]
            )
        )  # A/m2

        return np.sum(da * jboot, axis=0) / physics_variables.plascur

    @staticmethod
    def bootstrap_fraction_sakai(
        betap: float,
        q95: float,
        q0: float,
        alphan: float,
        alphat: float,
        eps: float,
        rli: float,
    ) -> float:
        """
        Calculate the bootstrap fraction using the Sakai formula.

        Parameters:
        betap (float): Plasma poloidal beta.
        q95 (float): Safety factor at 95% of the plasma radius.
        q0 (float): Safety factor at the magnetic axis.
        alphan (float): Density profile index
        alphat (float): Temperature profile index
        eps (float): Inverse aspect ratio.

        Returns:
        float: The calculated bootstrap fraction.

        Notes:
        The profile assumed for the alphan anf alpat indexes is only a prabolic profile without a pedestal (L-mode).
        The Root Mean Squared Error for the fitting database of this formula was 0.025
        Concentrating on the positive shear plasmas using the ACCOME code equilibria with the fully non-inductively driven
        conditions with neutral beam (NB) injection only are calculated.
        The electron temperature and the ion temperature were assumed to be equal
        This can be used for all apsect ratios.
        The diamagnetic fraction is included in this formula.

        References:
        Ryosuke Sakai, Takaaki Fujita, Atsushi Okamoto, Derivation of bootstrap current fraction scaling formula for 0-D system code analysis,
        Fusion Engineering and Design, Volume 149, 2019, 111322, ISSN 0920-3796,
        https://doi.org/10.1016/j.fusengdes.2019.111322.
        """
        # Sakai states that the ACCOME dataset used has the toridal diamagnetic current included in the bootstrap current
        # So the diamganetic current should not be calculated with this. idia = 0
        return (
            10 ** (0.951 * eps - 0.948)
            * betap ** (1.226 * eps + 1.584)
            * rli ** (-0.184 - 0.282)
            * (q95 / q0) ** (-0.042 * eps - 0.02)
            * alphan ** (0.13 * eps + 0.05)
            * alphat ** (0.502 * eps - 0.273)
        )

    def fhfac(self, is_):
        """Function to find H-factor for power balance
        author: P J Knight, CCFE, Culham Science Centre
        is : input integer : confinement time scaling law of interest
        This function calculates the H-factor required for power balance,
        using the given energy confinement scaling law.
        """

        physics_module.iscz = is_

        return root_scalar(self.fhz, bracket=(0.01, 150), xtol=0.003).root

    def fhz(self, hhh):
        """Function used to find power balance
        author: P J Knight, CCFE, Culham Science Centre
        hhh : input real : test value for confinement time H-factor
        This function is used to find power balance.
        <CODE>FHZ</CODE> is zero at power balance, which is achieved
        using routine <A HREF="zeroin.html">ZEROIN</A> to adjust the
        value of <CODE>hhh</CODE>, the confinement time H-factor.
        """

        (
            physics_variables.kappaa,
            ptrez,
            ptriz,
            taueezz,
            taueiz,
            taueffz,
            powerhtz,
        ) = self.pcond(
            physics_variables.afuel,
            physics_variables.palpmw,
            physics_variables.aspect,
            physics_variables.bt,
            physics_variables.dnitot,
            physics_variables.dene,
            physics_variables.dnla,
            physics_variables.eps,
            hhh,
            physics_variables.iinvqd,
            physics_module.iscz,
            physics_variables.ignite,
            physics_variables.kappa,
            physics_variables.kappa95,
            physics_variables.pchargemw,
            current_drive_variables.pinjmw,
            physics_variables.plascur,
            physics_variables.pcoreradpv,
            physics_variables.rmajor,
            physics_variables.rminor,
            physics_variables.te,
            physics_variables.ten,
            physics_variables.tin,
            physics_variables.q,
            physics_variables.qstar,
            physics_variables.vol,
            physics_variables.xarea,
            physics_variables.zeff,
        )

        # At power balance, fhz is zero.

        fhz = (
            ptrez
            + ptriz
            - physics_variables.falpha * physics_variables.palppv
            - physics_variables.pchargepv
            - physics_variables.pohmpv
        )

        # Take into account whether injected power is included in tau_e
        # calculation (i.e. whether device is ignited)

        if physics_variables.ignite == 0:
            fhz -= current_drive_variables.pinjmw / physics_variables.vol

        # Include the radiation power if requested

        if physics_variables.iradloss == 0:
            fhz += physics_variables.pradpv
        elif physics_variables.iradloss == 1:
            fhz += physics_variables.pcoreradpv

        return fhz

    @staticmethod
    def pcond(
        afuel,
        palpmw,
        aspect,
        bt,
        dnitot,
        dene,
        dnla,
        eps,
        hfact,
        iinvqd,
        isc,
        ignite,
        kappa,
        kappa95,
        pchargemw,
        pinjmw,
        plascur,
        pcoreradpv,
        rmajor,
        rminor,
        te,
        ten,
        tin,
        q,
        qstar,
        vol,
        xarea,
        zeff,
    ):
        """Routine to calculate the confinement times and
        the transport power loss terms.
        author: P J Knight, CCFE, Culham Science Centre
        afuel     : input real :  average mass of fuel (amu)
        palpmw    : input real :  alpha particle power (MW)
        aspect    : input real :  aspect ratio
        bt        : input real :  toroidal field on axis (T)
        dene      : input real :  volume averaged electron density (/m3)
        dnitot    : input real :  total ion density (/m3)
        dnla      : input real :  line-averaged electron density (/m3)
        eps       : input real :  inverse aspect ratio
        hfact     : input real :  H factor on energy confinement scalings
        iinvqd    : input integer :  switch for inverse quadrature
        isc       : input integer :  switch for energy confinement scaling to use
        ignite    : input integer :  switch for ignited calculation
        kappa     : input real :  plasma elongation
        kappa95   : input real :  plasma elongation at 95% surface
        kappaa    : output real : plasma elongation calculated using area ratio
        pchargemw : input real :  non-alpha charged particle fusion power (MW)
        pinjmw    : input real :  auxiliary power to ions and electrons (MW)
        plascur   : input real :  plasma current (A)
        pcoreradpv: input real :  total core radiation power (MW/m3)
        q         : input real :  edge safety factor (tokamaks), or
        rotational transform iotabar (stellarators)
        qstar     : input real :  equivalent cylindrical edge safety factor
        rmajor    : input real :  plasma major radius (m)
        rminor    : input real :  plasma minor radius (m)
        te        : input real :  average electron temperature (keV)
        ten       : input real :  density weighted average electron temp. (keV)
        tin       : input real :  density weighted average ion temperature (keV)
        vol       : input real :  plasma volume (m3)
        xarea     : input real :  plasma cross-sectional area (m2)
        zeff      : input real :  plasma effective charge
        ptrepv    : output real : electron transport power (MW/m3)
        ptripv    : output real : ion transport power (MW/m3)
        tauee     : output real : electron energy confinement time (s)
        taueff    : output real : global energy confinement time (s)
        tauei     : output real : ion energy confinement time (s)
        powerht   : output real : heating power (MW) assumed in calculation
        This subroutine calculates the energy confinement time
        using one of a large number of scaling laws, and the
        transport power loss terms.
        N. A. Uckan and ITER Physics Group,
        "ITER Physics Design Guidelines: 1989",
        ITER Documentation Series, No. 10, IAEA/ITER/DS/10 (1990)
        A. Murari et al 2015 Nucl. Fusion, 55, 073009
        C.C. Petty 2008 Phys. Plasmas, 15, 080501
        P.T. Lang et al. 2012 IAEA conference proceeding EX/P4-01
        ITER physics basis Chapter 2, 1999 Nuclear Fusion 39 2175
        Nuclear Fusion corrections, 2008 Nuclear Fusion 48 099801
        Menard 2019, Phil. Trans. R. Soc. A 377:20170440
        Kaye et al. 2006, Nucl. Fusion 46 848
        """

        eps2 = eps / 2.0e0
        str5 = 2.0e0 / (1.0e0 + (kappa**2))
        ck2 = (0.66e0 + (1.88e0 * (np.sqrt(eps2))) - (1.54e0 * eps2)) * (
            1.0e0 + (1.5e0 * (eps2**2))
        )
        chii = (
            (6.5e-22)
            * ck2
            * zeff
            * (aspect**1.5e0)
            * dene
            * (q**2)
            * str5
            / ((np.sqrt(tin)) * (bt**2))
        )
        str2 = 2.0e0 * (kappa**2) / (1.0e0 + (kappa**2))
        tauei = 0.375e0 * rminor**2 / chii * str2

        # Calculate heating power (MW)
        powerht = (
            physics_variables.falpha * palpmw + pchargemw + physics_variables.pohmmw
        )

        # If the device is not ignited, add the injected auxiliary power
        if ignite == 0:
            powerht = powerht + pinjmw

        # Include the radiation as a loss term if requested
        if physics_variables.iradloss == 0:
            powerht = powerht - physics_variables.pradpv * vol
        elif physics_variables.iradloss == 1:
            powerht = (
                powerht - pcoreradpv * vol
            )  # shouldn't this be vol_core instead of vol?
        # else do not adjust powerht for radiation

        # Ensure heating power is positive (shouldn't be necessary)
        powerht = max(powerht, 1.0e-3)

        # Line averaged electron density in scaled units
        dnla20 = dnla * 1.0e-20
        dnla19 = dnla * 1.0e-19

        # Volume averaged electron density in units of 10**20 m**-3
        n20 = dene / 1.0e20

        # Plasma current in MA
        pcur = plascur / 1.0e6

        # Separatrix kappa defined with X-section for general use
        kappaa = xarea / (np.pi * rminor * rminor)

        # Separatrix kappa defined with plasma volume for IPB scalings
        physics_variables.kappaa_ipb = vol / (
            2.0e0 * np.pi**2 * rminor * rminor * rmajor
        )

        # Calculate Neo-Alcator confinement time (used in several scalings)
        taueena = 0.07e0 * n20 * rminor * rmajor * rmajor * qstar

        # Electron energy confinement times

        if isc == 1:  # Neo-Alcator scaling (ohmic)
            # tauee = taueena
            tauee = hfact * taueena

        elif isc == 2:  # Mirnov scaling (H-mode)
            tauee = hfact * 0.2e0 * rminor * np.sqrt(kappa95) * pcur

        elif isc == 3:  # Merezhkin-Muhkovatov scaling (L-mode)
            tauee = (
                hfact
                * 3.5e-3
                * rmajor**2.75e0
                * rminor**0.25e0
                * kappa95**0.125e0
                * qstar
                * dnla20
                * np.sqrt(afuel)
                / np.sqrt(ten / 10.0e0)
            )

        elif isc == 4:  # Shimomura scaling (H-mode)
            tauee = (
                hfact
                * 0.045e0
                * rmajor
                * rminor
                * bt
                * np.sqrt(kappa95)
                * np.sqrt(afuel)
            )

        elif isc == 5:  # Kaye-Goldston scaling (L-mode)
            tauee = (
                hfact
                * 0.055e0
                * kappa95**0.28e0
                * pcur**1.24e0
                * n20**0.26e0
                * rmajor**1.65e0
                * np.sqrt(afuel / 1.5e0)
                / (bt**0.09e0 * rminor**0.49e0 * powerht**0.58e0)
            )

            if iinvqd != 0:
                tauee = 1.0e0 / np.sqrt(1.0e0 / taueena**2 + 1.0e0 / tauee**2)

        elif isc == 6:  # ITER Power scaling - ITER 89-P (L-mode)
            tauee = (
                hfact
                * 0.048e0
                * pcur**0.85e0
                * rmajor**1.2e0
                * rminor**0.3e0
                * np.sqrt(kappa)
                * dnla20**0.1e0
                * bt**0.2e0
                * np.sqrt(afuel)
                / np.sqrt(powerht)
            )

        elif isc == 7:  # ITER Offset linear scaling - ITER 89-O (L-mode)
            term1 = (
                0.04e0
                * pcur**0.5e0
                * rmajor**0.3e0
                * rminor**0.8e0
                * kappa**0.6e0
                * afuel**0.5e0
            )
            term2 = (
                0.064e0
                * pcur**0.8e0
                * rmajor**1.6e0
                * rminor**0.6e0
                * kappa**0.5e0
                * dnla20**0.6e0
                * bt**0.35e0
                * afuel**0.2e0
                / powerht
            )
            tauee = hfact * (term1 + term2)

        elif isc == 8:  # Rebut-Lallia offset linear scaling (L-mode)
            rll = (rminor**2 * rmajor * kappa95) ** 0.333e0
            tauee = (
                hfact
                * 1.65e0
                * np.sqrt(afuel / 2.0e0)
                * (
                    1.2e-2 * pcur * rll**1.5e0 / np.sqrt(zeff)
                    + 0.146e0
                    * dnla20**0.75e0
                    * np.sqrt(pcur)
                    * np.sqrt(bt)
                    * rll**2.75e0
                    * zeff**0.25e0
                    / powerht
                )
            )

        elif isc == 9:  # Goldston scaling (L-mode)
            tauee = (
                hfact
                * 0.037e0
                * pcur
                * rmajor**1.75e0
                * rminor ** (-0.37e0)
                * np.sqrt(kappa95)
                * np.sqrt(afuel / 1.5e0)
                / np.sqrt(powerht)
            )

            if iinvqd != 0:
                tauee = 1.0e0 / np.sqrt(1.0e0 / taueena**2 + 1.0e0 / tauee**2)

        elif isc == 10:  # T10 scaling
            denfac = dnla20 * rmajor * qstar / (1.3e0 * bt)
            denfac = min(1.0e0, denfac)
            tauee = (
                hfact
                * 0.095e0
                * rmajor
                * rminor
                * bt
                * np.sqrt(kappa95)
                * denfac
                / powerht**0.4e0
                * (
                    zeff**2
                    * pcur**4
                    / (rmajor * rminor * qstar**3 * kappa95**1.5e0)
                )
                ** 0.08e0
            )

        elif isc == 11:  # JAERI scaling
            gjaeri = (
                zeff**0.4e0
                * ((15.0e0 - zeff) / 20.0e0) ** 0.6e0
                * (
                    3.0e0
                    * qstar
                    * (qstar + 5.0e0)
                    / ((qstar + 2.0e0) * (qstar + 7.0e0))
                )
                ** 0.6e0
            )
            tauee = hfact * (
                0.085e0 * kappa95 * rminor**2 * np.sqrt(afuel)
                + 0.069e0
                * n20**0.6e0
                * pcur
                * bt**0.2e0
                * rminor**0.4e0
                * rmajor**1.6e0
                * np.sqrt(afuel)
                * gjaeri
                * kappa95**0.2e0
                / powerht
            )

        elif isc == 12:  # Kaye-Big scaling
            tauee = (
                hfact
                * 0.105e0
                * np.sqrt(rmajor)
                * rminor**0.8e0
                * bt**0.3e0
                * kappa95**0.25e0
                * pcur**0.85e0
                * n20**0.1e0
                * afuel**0.5e0
                / powerht**0.5e0
            )

        elif isc == 13:  # ITER H-mode scaling - ITER H90-P
            tauee = (
                hfact
                * 0.064e0
                * pcur**0.87e0
                * rmajor**1.82e0
                * rminor ** (-0.12e0)
                * kappa**0.35e0
                * dnla20**0.09e0
                * bt**0.15e0
                * np.sqrt(afuel)
                / np.sqrt(powerht)
            )

        elif isc == 14:  # Minimum of ITER 89-P (isc=6) and ITER 89-O (isc=7)
            tauit1 = (
                hfact
                * 0.048e0
                * pcur**0.85e0
                * rmajor**1.2e0
                * rminor**0.3e0
                * np.sqrt(kappa)
                * dnla20**0.1e0
                * bt**0.2e0
                * np.sqrt(afuel)
                / np.sqrt(powerht)
            )
            term1 = (
                0.04e0
                * pcur**0.5e0
                * rmajor**0.3e0
                * rminor**0.8e0
                * kappa**0.6e0
                * afuel**0.5e0
            )
            term2 = (
                0.064e0
                * pcur**0.8e0
                * rmajor**1.6e0
                * rminor**0.6e0
                * kappa**0.5e0
                * dnla20**0.6e0
                * bt**0.35e0
                * afuel**0.2e0
                / powerht
            )
            tauit2 = hfact * (term1 + term2)
            tauee = min(tauit1, tauit2)

        elif isc == 15:  # Riedel scaling (L-mode)
            tauee = (
                hfact
                * 0.044e0
                * pcur**0.93e0
                * rmajor**1.37e0
                * rminor ** (-0.049e0)
                * kappa95**0.588e0
                * dnla20**0.078e0
                * bt**0.152e0
                / powerht**0.537e0
            )

        elif isc == 16:  # Christiansen et al scaling (L-mode)
            tauee = (
                hfact
                * 0.24e0
                * pcur**0.79e0
                * rmajor**0.56e0
                * rminor**1.46e0
                * kappa95**0.73e0
                * dnla20**0.41e0
                * bt**0.29e0
                / (powerht**0.79e0 * afuel**0.02e0)
            )

        elif isc == 17:  # Lackner-Gottardi scaling (L-mode)
            qhat = (1.0e0 + kappa95**2) * rminor**2 * bt / (0.4e0 * pcur * rmajor)
            tauee = (
                hfact
                * 0.12e0
                * pcur**0.8e0
                * rmajor**1.8e0
                * rminor**0.4e0
                * kappa95
                * (1.0e0 + kappa95) ** (-0.8e0)
                * dnla20**0.6e0
                * qhat**0.4e0
                / powerht**0.6e0
            )

        elif isc == 18:  # Neo-Kaye scaling (L-mode)
            tauee = (
                hfact
                * 0.063e0
                * pcur**1.12e0
                * rmajor**1.3e0
                * rminor ** (-0.04e0)
                * kappa95**0.28e0
                * dnla20**0.14e0
                * bt**0.04e0
                * np.sqrt(afuel)
                / powerht**0.59e0
            )

        elif isc == 19:  # Riedel scaling (H-mode)
            tauee = (
                hfact
                * 0.1e0
                * np.sqrt(afuel)
                * pcur**0.884e0
                * rmajor**1.24e0
                * rminor ** (-0.23e0)
                * kappa95**0.317e0
                * bt**0.207e0
                * dnla20**0.105e0
                / powerht**0.486e0
            )

        elif isc == 20:  # Amended version of ITER H90-P law
            # Nuclear Fusion 32 (1992) 318
            tauee = (
                hfact
                * 0.082e0
                * pcur**1.02e0
                * bt**0.15e0
                * np.sqrt(afuel)
                * rmajor**1.60e0
                / (powerht**0.47e0 * kappa**0.19e0)
            )

        elif isc == 21:  # Large Helical Device scaling (stellarators)
            # S.Sudo, Y.Takeiri, H.Zushi et al., Nuclear Fusion 30 (1990) 11
            tauee = (
                hfact
                * 0.17e0
                * rmajor**0.75e0
                * rminor**2
                * dnla20**0.69e0
                * bt**0.84e0
                * powerht ** (-0.58e0)
            )

        elif isc == 22:  # Gyro-reduced Bohm scaling
            # R.J.Goldston, H.Biglari, G.W.Hammett et al., Bull.Am.Phys.Society,
            # volume 34, 1964 (1989)
            tauee = (
                hfact
                * 0.25e0
                * bt**0.8e0
                * dnla20**0.6e0
                * powerht ** (-0.6e0)
                * rminor**2.4e0
                * rmajor**0.6e0
            )

        elif isc == 23:  # Lackner-Gottardi stellarator scaling
            # K.Lackner and N.A.O.Gottardi, Nuclear Fusion, 30, p.767 (1990)
            iotabar = q  # dummy argument q is actual argument iotabar for stellarators
            tauee = (
                hfact
                * 0.17e0
                * rmajor
                * rminor**2
                * dnla20**0.6e0
                * bt**0.8e0
                * powerht ** (-0.6e0)
                * iotabar**0.4e0
            )

        elif (
            isc == 24
        ):  # ITER-93H scaling (ELM-free; multiply by 0.85 for ELMy version)
            # S.Kaye and the ITER Joint Central Team and Home Teams, in Plasma
            # Physics and Controlled Nuclear Fusion Research (Proc. 15th
            # Int. Conf., Seville, 1994) IAEA-CN-60/E-P-3
            tauee = (
                hfact
                * 0.053e0
                * pcur**1.06e0
                * bt**0.32e0
                * powerht ** (-0.67e0)
                * afuel**0.41e0
                * rmajor**1.79e0
                * dnla20**0.17e0
                * aspect**0.11e0
                * kappa**0.66e0
            )

        # Next two are ITER-97 H-mode scalings
        # J. G. Cordey et al., EPS Berchtesgaden, 1997

        elif isc == 26:  # ELM-free: ITERH-97P
            tauee = (
                hfact
                * 0.031e0
                * pcur**0.95e0
                * bt**0.25e0
                * powerht ** (-0.67e0)
                * dnla19**0.35e0
                * rmajor**1.92e0
                * aspect ** (-0.08e0)
                * kappa**0.63e0
                * afuel**0.42e0
            )

        elif isc == 27:  # ELMy: ITERH-97P(y)
            tauee = (
                hfact
                * 0.029e0
                * pcur**0.90e0
                * bt**0.20e0
                * powerht ** (-0.66e0)
                * dnla19**0.40e0
                * rmajor**2.03e0
                * aspect ** (-0.19e0)
                * kappa**0.92e0
                * afuel**0.2e0
            )

        elif isc == 28:  # ITER-96P (= ITER-97L) L-mode scaling
            # S.M.Kaye and the ITER Confinement Database Working Group,
            # Nuclear Fusion 37 (1997) 1303
            # N.B. tau_th formula used
            tauee = (
                hfact
                * 0.023e0
                * pcur**0.96e0
                * bt**0.03e0
                * kappa95**0.64e0
                * rmajor**1.83e0
                * aspect**0.06e0
                * dnla19**0.40e0
                * afuel**0.20e0
                * powerht ** (-0.73e0)
            )

        elif isc == 29:  # Valovic modified ELMy-H mode scaling
            tauee = (
                hfact
                * 0.067e0
                * pcur**0.9e0
                * bt**0.17e0
                * dnla19**0.45e0
                * afuel**0.05e0
                * rmajor**1.316e0
                * rminor**0.79e0
                * kappa**0.56e0
                * powerht ** (-0.68e0)
            )

        elif isc == 30:  # Kaye PPPL Workshop April 1998 L-mode scaling
            tauee = (
                hfact
                * 0.021e0
                * pcur**0.81e0
                * bt**0.14e0
                * kappa**0.7e0
                * rmajor**2.01e0
                * aspect ** (-0.18e0)
                * dnla19**0.47e0
                * afuel**0.25e0
                * powerht ** (-0.73e0)
            )

        elif isc == 31:  # ITERH-PB98P(y), ELMy H-mode scaling
            tauee = (
                hfact
                * 0.0615e0
                * pcur**0.9e0
                * bt**0.1e0
                * dnla19**0.4e0
                * powerht ** (-0.66e0)
                * rmajor**2
                * kappaa**0.75e0
                * aspect ** (-0.66e0)
                * afuel**0.2e0
            )

        elif isc == 32:  # IPB98(y), ELMy H-mode scaling
            # Data selection : full ITERH.DB3
            # Nuclear Fusion 39 (1999) 2175, Table 5
            tauee = (
                hfact
                * 0.0365e0
                * pcur**0.97e0
                * bt**0.08e0
                * dnla19**0.41e0
                * powerht ** (-0.63e0)
                * rmajor**1.93e0
                * kappa**0.67e0
                * aspect ** (-0.23e0)
                * afuel**0.2e0
            )

        elif isc == 33:  # IPB98(y,1), ELMy H-mode scaling
            # Data selection : full ITERH.DB3
            # Nuclear Fusion 39 (1999) 2175, Table 5
            tauee = (
                hfact
                * 0.0503e0
                * pcur**0.91e0
                * bt**0.15e0
                * dnla19**0.44e0
                * powerht ** (-0.65e0)
                * rmajor**2.05e0
                * physics_variables.kappaa_ipb**0.72e0
                * aspect ** (-0.57e0)
                * afuel**0.13e0
            )

        elif isc == 34:  # IPB98(y,2), ELMy H-mode scaling
            # Data selection : ITERH.DB3, NBI only
            # Nuclear Fusion 39 (1999) 2175, Table 5
            tauee = (
                hfact
                * 0.0562e0
                * pcur**0.93e0
                * bt**0.15e0
                * dnla19**0.41e0
                * powerht ** (-0.69e0)
                * rmajor**1.97e0
                * physics_variables.kappaa_ipb**0.78e0
                * aspect ** (-0.58e0)
                * afuel**0.19e0
            )

        elif isc == 35:  # IPB98(y,3), ELMy H-mode scaling
            # Data selection : ITERH.DB3, NBI only, no C-Mod
            # Nuclear Fusion 39 (1999) 2175, Table 5
            tauee = (
                hfact
                * 0.0564e0
                * pcur**0.88e0
                * bt**0.07e0
                * dnla19**0.40e0
                * powerht ** (-0.69e0)
                * rmajor**2.15e0
                * physics_variables.kappaa_ipb**0.78e0
                * aspect ** (-0.64e0)
                * afuel**0.20e0
            )

        elif isc == 36:  # IPB98(y,4), ELMy H-mode scaling
            # Data selection : ITERH.DB3, NBI only, ITER like devices
            # Nuclear Fusion 39 (1999) 2175, Table 5
            tauee = (
                hfact
                * 0.0587e0
                * pcur**0.85e0
                * bt**0.29e0
                * dnla19**0.39e0
                * powerht ** (-0.70e0)
                * rmajor**2.08e0
                * physics_variables.kappaa_ipb**0.76e0
                * aspect ** (-0.69e0)
                * afuel**0.17e0
            )

        elif isc == 37:  # ISS95 stellarator scaling
            # U. Stroth et al., Nuclear Fusion, 36, p.1063 (1996)
            # Assumes kappa = 1.0, triang = 0.0
            iotabar = q  # dummy argument q is actual argument iotabar for stellarators
            tauee = (
                hfact
                * 0.079e0
                * rminor**2.21e0
                * rmajor**0.65e0
                * dnla19**0.51e0
                * bt**0.83e0
                * powerht ** (-0.59e0)
                * iotabar**0.4e0
            )

        elif isc == 38:  # ISS04 stellarator scaling
            # H. Yamada et al., Nuclear Fusion, 45, p.1684 (2005)
            # Assumes kappa = 1.0, triang = 0.0
            iotabar = q  # dummy argument q is actual argument iotabar for stellarators
            tauee = (
                hfact
                * 0.134e0
                * rminor**2.28e0
                * rmajor**0.64e0
                * dnla19**0.54e0
                * bt**0.84e0
                * powerht ** (-0.61e0)
                * iotabar**0.41e0
            )

        elif isc == 39:  # DS03 beta-independent H-mode scaling
            # T. C. Luce, C. C. Petty and J. G. Cordey,
            # Plasma Phys. Control. Fusion 50 (2008) 043001, eqn.4.13, p.67
            tauee = (
                hfact
                * 0.028e0
                * pcur**0.83e0
                * bt**0.07e0
                * dnla19**0.49e0
                * powerht ** (-0.55e0)
                * rmajor**2.11e0
                * kappa95**0.75e0
                * aspect ** (-0.3e0)
                * afuel**0.14e0
            )

        elif isc == 40:  # "Non-power law" (NPL) Murari energy confinement scaling
            #  Based on the ITPA database of H-mode discharges
            #  A new approach to the formulation and validation of scaling expressions for plasma confinement in tokamaks
            #  A. Murari et al 2015 Nucl. Fusion 55 073009, doi:10.1088/0029-5515/55/7/073009
            #  Table 4.  (Issue #311)
            # Note that aspect ratio and M (afuel) do not appear, and B (bt) only
            # appears in the "saturation factor" h.
            h = dnla19**0.448e0 / (
                1.0e0 + np.exp(-9.403e0 * (bt / dnla19) ** 1.365e0)
            )
            tauee = (
                hfact
                * 0.0367e0
                * pcur**1.006e0
                * rmajor**1.731e0
                * kappaa**1.450e0
                * powerht ** (-0.735e0)
                * h
            )

        elif isc == 41:  # Beta independent dimensionless confinement scaling
            # C.C. Petty 2008 Phys. Plasmas 15, 080501, equation 36
            # Note that there is no dependence on the average fuel mass 'afuel'
            tauee = (
                hfact
                * 0.052e0
                * pcur**0.75e0
                * bt**0.3e0
                * dnla19**0.32e0
                * powerht ** (-0.47e0)
                * rmajor**2.09e0
                * kappaa**0.88e0
                * aspect ** (-0.84e0)
            )

        elif isc == 42:  # High density relevant confinement scaling
            # P.T. Lang et al. 2012, IAEA conference proceeding EX/P4-01
            # q should be q95: incorrect if icurr = 2 (ST current scaling)
            qratio = q / qstar
            # Greenwald density in m^-3
            nGW = 1.0e14 * plascur / (np.pi * rminor * rminor)
            nratio = dnla / nGW
            tauee = (
                hfact
                * 6.94e-7
                * plascur**1.3678e0
                * bt**0.12e0
                * dnla**0.032236e0
                * (powerht * 1.0e6) ** (-0.74e0)
                * rmajor**1.2345e0
                * physics_variables.kappaa_ipb**0.37e0
                * aspect**2.48205e0
                * afuel**0.2e0
                * qratio**0.77e0
                * aspect ** (-0.9e0 * np.log(aspect))
                * nratio ** (-0.22e0 * np.log(nratio))
            )

        elif isc == 43:  # Hubbard et al. 2017 I-mode confinement time scaling - nominal
            tauee = (
                hfact
                * 0.014e0
                * (plascur / 1.0e6) ** 0.68e0
                * bt**0.77e0
                * dnla20**0.02e0
                * powerht ** (-0.29e0)
            )

        elif isc == 44:  # Hubbard et al. 2017 I-mode confinement time scaling - lower
            tauee = (
                hfact
                * 0.014e0
                * (plascur / 1.0e6) ** 0.60e0
                * bt**0.70e0
                * dnla20 ** (-0.03e0)
                * powerht ** (-0.33e0)
            )

        elif isc == 45:  # Hubbard et al. 2017 I-mode confinement time scaling - upper
            tauee = (
                hfact
                * 0.014e0
                * (plascur / 1.0e6) ** 0.76e0
                * bt**0.84e0
                * dnla20**0.07
                * powerht ** (-0.25e0)
            )

        elif isc == 46:  # NSTX, ELMy H-mode scaling
            # NSTX scaling with IPB98(y,2) for other variables
            # Menard 2019, Phil. Trans. R. Soc. A 377:20170440
            # Kaye et al. 2006, Nucl. Fusion 46 848
            tauee = (
                hfact
                * 0.095e0
                * pcur**0.57e0
                * bt**1.08e0
                * dnla19**0.44e0
                * powerht ** (-0.73e0)
                * rmajor**1.97e0
                * physics_variables.kappaa_ipb**0.78e0
                * aspect ** (-0.58e0)
                * afuel**0.19e0
            )

        elif isc == 47:  # NSTX-Petty08 Hybrid
            # Linear interpolation between NSTX and Petty08 in eps
            # Menard 2019, Phil. Trans. R. Soc. A 377:20170440
            if (1.0e0 / aspect) <= 0.4e0:
                # Petty08, i.e. case (41)
                tauee = (
                    hfact
                    * 0.052e0
                    * pcur**0.75e0
                    * bt**0.3e0
                    * dnla19**0.32e0
                    * powerht ** (-0.47e0)
                    * rmajor**2.09e0
                    * kappaa**0.88e0
                    * aspect ** (-0.84e0)
                )

            elif (1.0e0 / aspect) >= 0.6e0:
                # NSTX, i.e.case (46)
                tauee = (
                    hfact
                    * 0.095e0
                    * pcur**0.57e0
                    * bt**1.08e0
                    * dnla19**0.44e0
                    * powerht ** (-0.73e0)
                    * rmajor**1.97e0
                    * physics_variables.kappaa_ipb**0.78e0
                    * aspect ** (-0.58e0)
                    * afuel**0.19e0
                )

            else:
                taupetty = (
                    0.052e0
                    * pcur**0.75e0
                    * bt**0.3e0
                    * dnla19**0.32e0
                    * powerht ** (-0.47e0)
                    * rmajor**2.09e0
                    * kappaa**0.88e0
                    * aspect ** (-0.84e0)
                )
                taunstx = (
                    0.095e0
                    * pcur**0.57e0
                    * bt**1.08e0
                    * dnla19**0.44e0
                    * powerht ** (-0.73e0)
                    * rmajor**1.97e0
                    * physics_variables.kappaa_ipb**0.78e0
                    * aspect ** (-0.58e0)
                    * afuel**0.19e0
                )

                tauee = hfact * (
                    (((1.0e0 / aspect) - 0.4e0) / (0.6e0 - 0.4e0)) * taunstx
                    + ((0.6e0 - (1.0e0 / aspect)) / (0.6e0 - 0.4e0)) * taupetty
                )

        elif isc == 48:  # NSTX gyro-Bohm (Buxton)
            # P F Buxton et al. 2019 Plasma Phys. Control. Fusion 61 035006
            tauee = (
                hfact
                * 0.21e0
                * pcur**0.54e0
                * bt**0.91e0
                * powerht ** (-0.38e0)
                * rmajor**2.14e0
                * dnla20 ** (-0.05e0)
            )

        elif isc == 49:  # tauee is an input
            tauee = hfact * physics_variables.tauee_in

        elif isc == 50:  # ITPA20 Issue #3164
            # The updated ITPA global H-mode confinement database: description and analysis
            # G. Verdoolaege et al 2021 Nucl. Fusion 61 076006 DOI 10.1088/1741-4326/abdb91

            # thermal energy confinement time
            # plasma current Ip (MA),
            # on-axis vacuum toroidal magnetic field Bt (T)
            # "central line-averaged electron density" nebar (1019 m−3)
            # thermal power lost due to transport through the LCFS Pl,th (MW)
            # major radius Rgeo (m)
            # elongation of the LCFS, defined as κa = V/(2πRgeo πa2)
            # (with V (m3) the plasma volume inside the LCFS),
            # inverse aspect ratio epsilon = a/Rgeo
            # effective atomic mass Meff of the plasma - NOT defined, but I have taken it equal to
            # aion = average mass of all ions (amu).
            # energy confinement time is given by τE,th = Wth/Pl,th, where Wth is the thermal stored energy.
            # The latter is derived from the total stored energy Wtot by subtracting the energy
            # content associated to fast particles originating from plasma heating.

            tauee = (
                hfact
                * 0.053
                * pcur**0.98
                * bt**0.22
                * dnla19**0.24
                * powerht ** (-0.669)
                * rmajor**1.71
                * (1 + physics_variables.triang) ** 0.36
                * physics_variables.kappaa_ipb**0.8
                * eps**0.35
                * physics_variables.aion**0.2
            )

        else:
            error_handling.idiags[0] = isc
            error_handling.report_error(81)

        # Ion energy confinement time
        # N.B. Overwrites earlier calculation above

        tauei = tauee

        # Calculation of the transport power loss terms
        # Transport losses in Watts/m3 are 3/2 * n.e.T / tau , with T in eV
        # (here, tin and ten are in keV, and ptrepv and ptripv are in MW/m3)

        ptripv = 2.403e-22 * dnitot * tin / tauei
        ptrepv = 2.403e-22 * dene * ten / tauee

        ratio = dnitot / dene * tin / ten

        # Global energy confinement time

        taueff = (ratio + 1.0e0) / (ratio / tauei + 1.0e0 / tauee)

        return kappaa, ptrepv, ptripv, tauee, tauei, taueff, powerht
