"""
Library functions for the density and temperature profiles in PROCESS
Also contain residual functions for fitting

Author: Hanni Lux (Hanni.Lux@ccfe.ac.uk)

Compatible with PROCESS version 379
"""


from scipy.special import gamma
from numpy import array, diff, pi, sin, argmin

NUMACC = 1e-10


def ncore(rhopedn, nped, nsep, nav, alphan):

    """Calculates the core density from the average density
    and other profile parameters"""

    return (
        (1.0 / 3.0)
        * (1.0 / rhopedn**2)
        * (
            3.0 * nav * (1.0 + alphan)
            + nsep * (1.0 + alphan) * (-2.0 + rhopedn + rhopedn**2)
            - nped
            * (
                1.0
                + alphan
                + rhopedn
                + (alphan * rhopedn)
                + (alphan - 2.0) * rhopedn**2
            )
        )
    )


def nprofile(rho, rhopedn, nzero, nped, nsep, alphan):

    """Calculates the density at a given normalised radius rho
    depending on the profile parameters"""

    if rho < 0.0:
        print("Error: rho undefined!")
    elif rho <= rhopedn:
        return nped + (nzero - nped) * (1.0 - (rho**2 / rhopedn**2)) ** alphan
    elif rho <= 1.0:
        return nsep + (nped - nsep) * (1.0 - rho) / (1.0 - rhopedn)
    else:
        print("Error: rho undefined!")


def tcore(rhopedt, tped, tsep, tav, alphat, tbeta):

    """Calculates the core temperature from the average temperature
    and other profile parameters"""

    gamfac = (
        gamma(1.0e0 + alphat + 2.0e0 / tbeta)
        / gamma((2.0e0 + tbeta) / tbeta)
        / rhopedt**2.0e0
    )

    if (abs(alphat % 1) <= NUMACC) or (abs((alphat % 1) - 1) <= NUMACC):
        gamfac = -gamfac / gamma(1.0e0 + alphat)
    else:
        gamfac = gamfac * gamma(-alphat) * sin(pi * alphat) / pi

    return (
        tped
        + (
            tped * rhopedt**2.0e0
            - tav
            + (1.0e0 - rhopedt)
            / 3.0e0
            * ((1.0e0 + 2.0e0 * rhopedt) * tped + (2.0e0 + rhopedt) * tsep)
        )
        * gamfac
    )


def tprofile(rho, rhopedt, tzero, tped, tsep, alphat, betat):

    """Calculates the temperature at a given normalised radius rho
    depending on the profile parameters"""

    if rho < 0.0:
        print("Error: rho undefined!")
    elif rho <= rhopedt:
        return (
            tped + (tzero - tped) * (1.0 - (rho**betat / rhopedt**betat)) ** alphat
        )
    elif rho <= 1.0:
        return tsep + (tped - tsep) * (1.0 - rho) / (1.0 - rhopedt)
    else:
        print("Error: rho undefined!")


def nresidual(params, xdata, ydata, dictfitrad):

    """Calculates the residual between a given data set and a given
    density profile"""

    alphan, rhoped = params

    nzero = dictfitrad["n0"]
    indped = argmin(abs(xdata - rhoped))
    dictfitrad["nped"] = ydata[indped]
    nped = dictfitrad["nped"]
    nsep = dictfitrad["nsep"]

    ymodel = [nprofile(x, rhoped, nzero, nped, nsep, alphan) for x in xdata]

    return array(ymodel) - ydata


def tresidual(params, xdata, ydata, dictfitrad):

    """Calculates the residual between a given data set and a given
    temperature profile"""

    alphat, betat, rhoped = params

    tzero = dictfitrad["t0"]
    indped = argmin(abs(xdata - rhoped))
    dictfitrad["tped"] = ydata[indped]
    tped = dictfitrad["tped"]
    tsep = dictfitrad["tsep"]

    ymodel = [tprofile(x, rhoped, tzero, tped, tsep, alphat, betat) for x in xdata]

    return array(ymodel) - ydata


def int_f_rho_drho(rhoarr, parr):

    """integrates f * rho drho using a leftsided rectangular sum"""

    ptot = 0.0

    drhoarr = diff(rhoarr)

    for rho, p, drho in zip(rhoarr[:-1], parr[:-1], drhoarr):

        ptot += p * rho * drho

    return ptot


if __name__ == "__main__":

    RHOPEDN = 0.85
    NPED = 1e2
    ALPHAN = 2.0
    NAV = 1.4e2
    NSEP = 0.42
    NCORE = ncore(RHOPEDN, NPED, NSEP, NAV, ALPHAN)

    from pylab import figure, xlabel, ylabel, axvline, axhline, plot, show
    from numpy import linspace

    RHO = linspace(0, 1, 50)

    figure()
    xlabel(r"$\rho$")
    ylabel(r"$n$")
    axvline(RHOPEDN, label=r"$\rho_{ped}$")
    axhline(NCORE, ls="--")
    axhline(NPED, ls=":")
    axhline(NSEP)
    plot(RHO, [nprofile(x, RHOPEDN, NCORE, NPED, NSEP, ALPHAN) for x in RHO])

    show()
