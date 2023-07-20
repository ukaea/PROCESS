#!/usr/bin/env python
"""

POPCON (Plasma OPeration CONtour) Plot
       (Power OutPut CONtour) Plot

Stuart Muldrew (stuart.muldrew@ukaea.uk)
Updated version of code by PJ Knight and TN Todd
12/09/2018

Outputs POPCON plot for given input parameters.

For reference, density and temperature profiles are of the form:

n(r) = n0 (1 - (r/a)**2)**alphan   (similarly temperature with index alphat)

Line averaged value:    nbar = n0 / sqrt(1.0 + alphan)

Volume averaged value:   <n> = n0 / (1.0 + alphan)
                   or:   <n> = nbar / sqrt(1.0 + alphan)


History:
- Originally written by TNT in late 1980s (BBC Microcomputer / BBC Basic / primitive graphics)
- Ported to IDL by PJK in 2011, added new routines and compared with W Han's cut-down PROCESS
- Ported to Python in 2018 by SIM and added MFILE reading


"""
# Imported libraries
import argparse
import process.io.mfile as mf
import numpy as np
import matplotlib.pyplot as plt
import sys


# Functions
def ohmic():
    """Ohmic Power"""
    Foh = (
        plascur
        * plascur
        / 16.0
        * R
        / (a * a * K)
        * Zeff
        * (1.0 + 1.5 * alphat)
        / fnc
        / Trat**1.5
    )
    Poh = Foh / T0**1.5

    return Poh


def brem_tnt():
    """Bremsstrahlung - TNT model"""
    fPb = Zeff * np.sqrt(Trat) * fv / (1.0 + 2.0 * alphan + 0.5 * alphat) / 1000.0
    Pbf = n0 * n0 * fPb
    Pb = Pbf * np.sqrt(T0)

    return Pb


def brem_process():
    """Bremsstrahlung - Process Model"""
    kappa95 = 0.96 * K  # RUN66 ratio of kappa95 to kappaa
    vol = 20.0 * fv

    dene = 1.0e19 * nbar / np.sqrt(1.0 + alphan)
    den20 = dene / 1.0e20

    pcoef = (1.0 + alphan) * (1.0 + alphat) / (1.0 + alphan + alphat)
    ten = T0 / (1.0 + alphat) * pcoef
    t10 = ten / 10.0

    fbhe = 0.9
    impc = 1.0
    fbc = 0.52
    impo = 1.0
    fbo = 0.52
    impfe = 1.0
    fbfe = 0.35
    cfe0 = 0.35415e-3

    ralpne = 0.1
    rncne = impc * (0.009 + 0.006 * (7.0e19 / dene) ** 2.6)
    rnone = impo * 0.001
    rnfene = impfe * (0.0005 * (7.0e19 / dene) ** 2.3 + cfe0)

    radexp = (
        (1.0 + alphan) ** 1.5
        * np.sqrt(1.0 + alphan + alphat)
        / (1.0 + 2.0 * alphan + 0.5 * alphat)
    )

    deni = dene * (
        1.0 - 2.0 * ralpne - (6.0 * rncne + 8.0 * rnone + 26.0 * rnfene)
    )  # No beam

    vr = R * (a * (1.0 + kappa95) / 2.0) ** 2 / (58.652 * vol)
    phe = 65.8 * ralpne * (dene / 7.0e19) ** 1.5 * vr
    pc = 1120.0 * rncne * (dene / 7.0e19) ** 1.5 * vr
    po = 2240.0 * rnone * (dene / 7.0e19) ** 1.5 * vr
    pfe = 44800.0 * rnfene * (dene / 7.0e19) ** 2.5 * vr

    pbremdt = 0.016 * radexp * den20 * den20 * (deni / dene) * np.sqrt(t10)  # D-T Brem
    pbremz = fbhe * phe + fbc * pc + fbo * po + fbfe * pfe  # High Z Bremsstrahlung

    return (pbremz + pbremdt) * vol


def sync_tnt():
    """Synchrotron - TNT model"""
    Psf1 = 18.0 * a / R / np.sqrt(Trat)
    fPs = (
        20.0
        * fv
        * 4.1e-7
        * (Trat * B) ** 2.5
        * np.sqrt((1.0 - refl) / 10.0 / a / np.sqrt(K))
    )
    Psf = fPs * np.sqrt(n0)
    Psf5 = fion / 250.0 / (B * B) / (1.0 + 2.0 * alphan + 2.0 * alphat)
    s1 = Psf * T0**2.5 * np.sqrt(1.0 + Psf1 / np.sqrt(T0))
    s5 = n0 * T0 * Psf5
    Ps = s1 * (1.0 / fp + s4 - s5)

    return Ps


def sync_process():
    """Synchrotron - Process Model"""
    pcoef = (1.0 + alphan) * (1.0 + alphat) / (1.0 + alphan + alphat)
    ten = T0 / (1.0 + alphat) * pcoef / 10.0
    den20 = nbar / np.sqrt(1.0 + alphan) / 10.0
    xfact = 5.7 * a / R / np.sqrt(ten)

    Ps = (
        1.3e-4
        * (B * ten) ** 2.5
        * np.sqrt(1.0 + xfact)
        * np.sqrt(den20 / a)
        * np.sqrt(1.0 - refl)
        * fv
        * 20.0
    )
    # MW/m3 to MW conversion

    return Ps


def sync_albajar_fidone():
    """Synchroton - Albajar and Fidone Model"""
    tbet = 2.0  # tbet is betaT in Albajar, not to be confused with plasma beta
    rpow = 0.62  # rpow is the (1-Rsyn) power dependence based on plasma shape

    de2o = n0 / 10.0
    teo = T0 * Trat
    pao = 6040.0 * (a * de2o) / B
    gfun = 0.93 * (1.0 + 0.85 * np.exp(-0.82 * R / a))
    kfun = (alphan + 3.87 * alphat + 1.46) ** (-0.79)
    kfun = kfun * (1.98 + alphat) ** 1.36 * tbet**2.14
    kfun = kfun * (tbet**1.53 + 1.87 * alphat - 0.16) ** (-1.33)
    dum = 1.0 + 0.12 * (teo / (pao**0.41)) * (1.0 - refl) ** 0.41

    dum = dum ** (-1.51)  # Very high T modification, from Fidone

    Ps = 3.84e-8 * (1.0 - refl) ** rpow * R * a**1.38
    Ps = Ps * K**0.79 * B**2.62 * de2o**0.38
    Ps = Ps * teo * (16.0 + teo) ** 2.61 * dum * gfun * kfun

    return Ps


def fusion():
    """Charged particle fusion power"""
    if fuel2 == "T":
        FPalf = 1.52e-6 * fv * n0 * n0 * ((Zimp - Zeff) / (Zimp - 1.0)) ** 2
        g = 23.9 - 17.87 * T0**0.04
        GGG = 23.9 - (17.37 + 0.7 * np.log(T0)) * T0**0.04
    else:
        FPalf = 1.16e-7 * fv * n0 * n0 * ((Zimp - Zeff) / (1.8 * Zimp - 2.6)) ** 2
        g = 14.81 - 10.0 * T0**0.03
        GGG = 14.81 - (9.79 + 0.3 * np.log(T0)) * T0**0.03

    Palf = FPalf * (T0**g) / (1.0 + 2.0 * alphan + GGG * alphat)

    return Palf


def transport():
    """
    Calculate transport power given the stored plasma thermal energy
    and a chosen scaling law for the confinement time

    etherm = plasma thermal energy in MJ
    nbar = line-averaged density in 1.0E19 /m3
    afuel = average fuel mass (amu)

    ptransp = etherm/taue, but taue = ftaue*(Ptransp**pheat_index),
    hence the coding below
    """

    if scaling_law == 0:  # IPB98(y,2) scaling
        etherm = 0.4 * C1 / (B * B) * n0 / fp * fion * T0

        # ftaue is the confinement time scaling, but excluding the pheat term
        ftaue = (
            hfact
            * 0.0562
            * plascur**0.93
            * B**0.15
            * nbar**0.41
            * R**1.39
            * a ** (0.58)
            * K**0.78
            * afuel**0.19
        )

        pheat_index = -0.69

        Ptransp = (etherm / ftaue) ** (1.0 / (1.0 + pheat_index))

    else:
        print("scaling_law =", scaling_law)
        print("No such scaling law!")
        sys.exit()

    return Ptransp


if __name__ == "__main__":

    # Setup command line arguments
    parser = argparse.ArgumentParser(
        description="Displays a POPCON plot for input MFILE.  "
        "For more information contact Stuart.Muldrew@ukaea.uk or Peter.Knight@ukaea.uk"
    )

    parser.add_argument(
        "-f",
        metavar="MFILE",
        type=str,
        default="MFILE.DAT",
        help="specify the MFILE (default=MFILE.DAT)",
    )

    parser.add_argument(
        "-s", "--save", help="save as well as displaying figure", action="store_true"
    )

    parser.add_argument(
        "-t",
        "--test",
        help="Test mode: ignores MFILE and runs with R Kembleton's test values",
        action="store_true",
    )

    parser.add_argument(
        "-x",
        type=int,
        default="1",
        help="Temperature (x-axis): (0) Central (1) Volume (default=1)",
    )

    parser.add_argument(
        "-y",
        type=int,
        default="1",
        help="Density (y-axis): (0) Central (1) Volume (2) Line (default=1)",
    )

    parser.add_argument(
        "-z",
        type=int,
        default="0",
        help="Power (z-axis): (0) Aux for balance (1) Net (2) Fusion (default=0)",
    )

    parser.add_argument(
        "-zimp", type=float, default="18.0", help="Impurity charge (default=18 Ar)"
    )

    args = parser.parse_args()

    # Option for Test mode that ignores MFILE and runs with set parameters
    if args.test:
        # Test parameters from Richard Kembleton's DEMO1 (01/10/2013)
        a = 2.2  # Minor radius (m)
        R = 9.0  # Major radius (m)
        K = 1.8  # Elongation
        plascur = 12.0  # Plasma current (MA)
        refl = 0.6  # Wall reflection coeff.
        B = 7.33  # Toroidal magnetic field
        alphan = 0.37  # Density profile index
        alphat = 1.63  # Temperature profile index
        Zimp = 18.0  # Impurity charge
        Zeff = 3.48  # Zeff
        Trat = 1.0  # Ratio of Te/Ti
        Paux = 70.0  # Auxiliary power (MW)
        fuel2 = "T"  # T or He3 species mix with D
        hfact = 1.4  # Confinement time H factor
        j0fudge = 0.367  # Only used for Ohmic power calculation
        RadModel = 2  # Radiation Model: (0) TNT (1) Process (2) Process + Albajar and Fidone Sync Model
    else:
        m_file = mf.MFile(args.f)
        a = m_file.data["rminor"].get_scan(-1)  # Minor radius (m)
        R = m_file.data["rmajor"].get_scan(-1)  # Major radius (m)
        K = m_file.data["kappaa"].get_scan(-1)  # Elongation
        plascur = m_file.data["plascur/1d6"].get_scan(-1)  # Plasma current (MA)
        refl = m_file.data["ssync"].get_scan(-1)  # Wall reflection coeff.
        B = m_file.data["bt"].get_scan(-1)  # Toroidal magnetic field
        alphan = m_file.data["alphan"].get_scan(-1)  # Density profile index
        alphat = m_file.data["alphat"].get_scan(-1)  # Temperature profile index
        Zeff = m_file.data["zeff"].get_scan(-1)  # Zeff
        te = m_file.data["te"].get_scan(-1)  # Electron temperature
        ti = m_file.data["ti"].get_scan(-1)  # Ion temperature
        Paux = m_file.data["pinjmw"].get_scan(-1)  # Auxiliary power (MW)
        hfact = m_file.data["hfact"].get_scan(-1)  # Confinement time H factor

        Zimp = args.zimp  # Impurity charge

        Trat = te / ti  # Ratio of Te/Ti
        fuel2 = "T"  # T or He3 species mix with D
        j0fudge = 0.367  # Only used for Ohmic power calculation
        RadModel = 2  # Radiation Model: (0) TNT (1) Process (2) Process + Albajar and Fidone Sync Model

    # Array Size
    nmin = 2.0976
    nmax = 25.0
    tmin = 5.0
    tmax = 40.0

    tpoints = 120
    npoints = 100

    # Plot
    output = "popcon.pdf"
    # Temperature (x-axis): (0) Central (1) Volume
    temp = args.x
    # Density (y-axis): (0) Central (1) Volume (2) Line
    dens = args.y
    # Power (z-axis): (0) Aux for balance (1) Net
    Pplot = args.z

    # Scaling (Only one option)
    scaling_law = 0
    scaling_name = "IPB98(y,2)"

    #####################################

    # Main Program

    # Imports

    if fuel2 == "He3":
        if Zeff == 1.445:
            Zeff = 1.445

    # Critical density, I/N scaling

    ncrit = 6.67 * plascur / (np.pi * a * a * K)

    if nmax == 0.0:
        nmax = ncrit

    # Line-averaged electron density, 1.0E19 m-3
    nbar = nmin + np.arange(npoints) / (npoints - 1) * (nmax - nmin)
    nbar = np.flip(nbar, 0)

    # Central ion temperature, keV
    T0 = tmin + np.arange(tpoints) / (tpoints - 1) * (tmax - tmin)

    # Meshgrid
    T0, nbar = np.meshgrid(T0, nbar)

    # Perform main calculation at the given line-averaged density, central temperature
    n0 = nbar * np.sqrt(1.0 + alphan)  # Central density

    fnc = (1.0 - np.sqrt(2.63 * a / R / (5.0 + alphat))) ** 2
    j0 = j0fudge * plascur / (np.pi * a * a * K) * (1.0 + 1.5 * alphat) / fnc

    fv = a * a * R * K  # Volume factor
    C1 = B * B * fv / 8.49
    fp = 1.0 + alphan + alphat  # Profile factor
    s3 = fp + 0.72 * (1.0 + 3.0 * alphat)
    s4 = (
        4.16
        * ((1.0 + 1.5 * alphat) * plascur / 5.0 / a / B / K) ** 2
        / (2.0 + 4.5 * alphat)
        / s3
    )

    if fuel2 == "T":
        fion = Trat + (1.0 + Zimp - Zeff) / Zimp
        nrat = (Zeff - 1.0) / Zimp / (Zimp - Zeff)
        afuel = (2.5 + 2.0 * nrat) / (1.0 + nrat)
    else:
        fion = Trat + (Zimp * (Zimp - Zeff) + 1.8 * Zeff - 2.6) / Zimp / (
            1.8 * Zimp - 2.6
        )
        nrat = (1.8 * Zeff - 2.6) / Zimp / (Zimp - Zeff)
        afuel = 2.0 - 0.2 / (1.0 + nrat)

    # Power Calculations
    Pohm = ohmic()

    if RadModel == 0:
        Pbrem = brem_tnt()
        Psync = sync_tnt()
    elif RadModel == 1:
        Pbrem = brem_process()
        Psync = sync_process()
    elif RadModel == 2:
        Pbrem = brem_process()
        Psync = sync_albajar_fidone()

    Palf = fusion()
    Ptransp = transport()

    # Total (continuum) radiation power
    Prad = Pbrem + Psync
    # Total heating power (excluding Paux)
    Pheat = Palf + Pohm

    # Auxiliary power required for power balance (MW)
    Pbal = -(Pheat - Ptransp - Prad)

    # Net power output (MW)
    Pnet = Pheat + Paux - Ptransp - Prad

    #
    if temp == 0:
        xtitle = "Central ion temperature / keV"
        xplot = T0[0, :]
    elif temp == 1:
        xtitle = "Volume averaged ion temperature / keV"
        xplot = T0[0, :] / (1.0 + alphat)

    if dens == 0:
        ytitle = r"Central electron density / $10^{19} \mathrm{m}^3$"
        yplot = nbar[:, 0] * np.sqrt(1.0 + alphan)
    elif dens == 1:
        ytitle = r"Volume averaged electron density / $10^{19} \mathrm{m}^3$"
        yplot = nbar[:, 0] / np.sqrt(1.0 + alphan)
    elif dens == 2:
        ytitle = r"Line-averaged electron density / $10^{19} \mathrm{m}^3$"
        yplot = nbar[:, 0]

    # Plotting
    if Pplot == 0:
        plt.contourf(
            Pbal,
            60,
            extent=[np.min(xplot), np.max(xplot), np.max(yplot), np.min(yplot)],
            cmap="jet",
        )
        plt.colorbar(label="$P_{aux}$ / MW")
        plt.contour(
            Pbal,
            levels=[0.0],
            extent=[np.min(xplot), np.max(xplot), np.max(yplot), np.min(yplot)],
            colors="k",
        )
    elif Pplot == 1:
        plt.contourf(
            Pnet,
            60,
            extent=[np.min(xplot), np.max(xplot), np.max(yplot), np.min(yplot)],
            cmap="jet",
        )
        plt.colorbar(label="$P_{net}$ / MW")
        plt.contour(
            Pnet,
            levels=[0.0],
            extent=[np.min(xplot), np.max(xplot), np.max(yplot), np.min(yplot)],
            colors="k",
        )
    elif Pplot == 2:
        plt.contourf(
            5.0 * Palf,
            60,
            extent=[np.min(xplot), np.max(xplot), np.max(yplot), np.min(yplot)],
            cmap="jet",
        )
        plt.colorbar(label="$P_{fus}$ / MW")

    plt.xlabel(xtitle)
    plt.ylabel(ytitle)
    if args.save:
        plt.savefig(output)
    plt.show()
