#!/usr/bin/env python
"""
Python Tool to fit a general temperature or density profile from an ascii table
using our pedestal parametrisation

Author: H. Lux (Hanni.Lux@ccfe.ac.uk)

Input file:
profile.txt (can be either density or temperature profile)
(should contain data in columns, comment lines start with #)

Compatible with PROCESS version 379
"""

#######################
# imported libraries

import argparse
from numpy import loadtxt, array, argsort, arange, argmin, linspace
from scipy.optimize import leastsq
from process_io_lib.profile_funcs import nresidual, tresidual, nprofile, tprofile
from pylab import figure, xlabel, ylabel, axvline, axhline, plot, show

if __name__ == "__main__":
    ############################################################
    # Usage

    PARSER = argparse.ArgumentParser(
        description="Program to fit a\
 general temperature or density profile from an ascii table."
    )

    PARSER.add_argument(
        "-f",
        "--filename",
        default="profile.txt",
        help="ascii file containing data in columns,\
 default = profile.txt",
    )

    PARSER.add_argument(
        "-r",
        "--rho",
        type=int,
        default=0,
        help="column of the normalised radius rho=r/a,\
 default = 0",
    )

    PARSER.add_argument(
        "-n",
        "--density",
        type=int,
        default=1,
        help="column of the density profile,\
 default = 1",
    )

    PARSER.add_argument(
        "-t",
        "--temperature",
        type=int,
        default=2,
        help="column of the temperature profile,\
 default = 2",
    )

    PARSER.add_argument(
        "-rn",
        "--rhopedn",
        type=float,
        default=0.9,
        help="user defined initial guess of the density\
 pedestal position, if outside [0,1] starts at 0.9, default = 0.9",
    )

    PARSER.add_argument(
        "-rt",
        "--rhopedt",
        type=float,
        default=0.9,
        help="user defined initial guess of the temperature\
 pedestal position, if outside [0,1], starts at 0.9, default = 0.9",
    )

    ARGS = PARSER.parse_args()

    USEDENSITY = True
    USETEMPERATURE = True

    #####################################
    # input error handling

    try:
        DATA = loadtxt(ARGS.filename, unpack=True)
    except FileNotFoundError:
        print("Error: There is no file called", ARGS.filename)
        exit()
    except ValueError:
        print(
            "Error: In",
            ARGS.filename,
            "all comment rows must\
 start with a #! Columns can only contain floats!",
        )
        exit()

    if DATA.size == 0:
        print("Error: There is no data in", ARGS.filename)
        exit()

    try:
        RHO = DATA[ARGS.rho]
    except IndexError:
        print(
            "Error: The column for the normalised radius does not exist!\
 Remember to start counting at 0!"
        )
        exit()

    try:
        DEN = DATA[ARGS.density]
    except IndexError:
        print(
            "Warning: The column for the density does not exist!\
 Remember to start counting at 0!"
        )
        USEDENSITY = False

    try:
        TEMP = DATA[ARGS.temperature]
    except IndexError:
        print(
            "Warning: The column for the temperature does not exist!\
 Remember to start counting at 0!"
        )
        USETEMPERATURE = False

    ###################################################
    # assure sorted data
    SORTED_INDS = argsort(RHO)

    if all(SORTED_INDS != arange(len(RHO))):
        RHO = RHO[SORTED_INDS]

        if USEDENSITY:
            DEN = DEN[SORTED_INDS]

        if USETEMPERATURE:
            TEMP = TEMP[SORTED_INDS]

    # only use data between 0 and 1!
    # TODO: Does not work, if values outside range are closer!
    # However, this is better, for getting most accurate values at rho = 0,1!
    INDZERO = argmin(abs(RHO))
    INDONE = argmin(abs(RHO - 1.0))
    RHO = RHO[INDZERO : INDONE + 1]

    if USEDENSITY:
        DEN = DEN[INDZERO : INDONE + 1]

    if USETEMPERATURE:
        TEMP = TEMP[INDZERO : INDONE + 1]

    ###################################################

    if USEDENSITY:

        DICTFIT = {"alphan": 2.0}

        if ARGS.rhopedn >= 0.0 and ARGS.rhopedn <= 1.0:
            DICTFIT["rhopedn"] = ARGS.rhopedn
        else:
            DICTFIT["rhopedn"] = 0.9
        DICTFIT["n0"] = DEN[0]
        DICTFIT["nsep"] = DEN[-1]
        INDPED = argmin(abs(RHO - DICTFIT["rhopedn"]))
        DICTFIT["nped"] = DEN[INDPED]

        LSQFIT = leastsq(
            nresidual,
            array([DICTFIT["alphan"], DICTFIT["rhopedn"]]),
            args=(RHO, DEN, DICTFIT),
        )
        if LSQFIT[-1] not in [1, 2, 3, 4]:
            print("Error: Density fit has failed!")
            print("leastsq return:", LSQFIT)

        else:
            print("\nDensity profile parameters:")
            print("----------------------------")
            print("Density pedestal rhopedn:", LSQFIT[0][1])
            print("Core density n0:", DICTFIT["n0"])
            print("Pedestal density nped:", DICTFIT["nped"])
            print("Density at separatrix nsep:", DICTFIT["nsep"])
            print("Density peaking parameter alphan:", LSQFIT[0][0], "\n")

            RHOARR = linspace(0, 1, 100)
            RHOPEDN = LSQFIT[0][1]
            NCORE = DICTFIT["n0"]
            NPED = DICTFIT["nped"]
            NSEP = DICTFIT["nsep"]
            ALPHAN = LSQFIT[0][0]

            figure()
            xlabel(r"$\rho$")
            ylabel(r"$n$")
            axvline(RHOPEDN, label=r"$\rho_{ped}$")
            axhline(NCORE, ls="--")
            axhline(NPED, ls=":")
            axhline(NSEP)
            plot(
                RHOARR,
                [nprofile(x, RHOPEDN, NCORE, NPED, NSEP, ALPHAN) for x in RHOARR],
            )
            plot(RHO, DEN, "kx")

    if USETEMPERATURE:
        DICTFIT = {"alphat": 2.0, "betat": 2.0}

        if ARGS.rhopedt >= 0.0 and ARGS.rhopedt <= 1.0:
            DICTFIT["rhopedt"] = ARGS.rhopedt
        else:
            DICTFIT["rhopedt"] = 0.9
        DICTFIT["t0"] = TEMP[0]
        DICTFIT["tsep"] = TEMP[-1]
        INDPED = argmin(abs(RHO - DICTFIT["rhopedt"]))
        DICTFIT["tped"] = TEMP[INDPED]

        LSQFIT = leastsq(
            tresidual,
            array([DICTFIT["alphat"], DICTFIT["betat"], DICTFIT["rhopedt"]]),
            args=(RHO, TEMP, DICTFIT),
        )

        if LSQFIT[-1] not in [1, 2, 3, 4]:
            print("Error: Temperature fit has failed!")
            print("leastsq return:", LSQFIT)

        else:
            print("\nTemperature profile parameters:")
            print("----------------------------")
            print("Temperature pedestal rhopedt:", LSQFIT[0][2])
            print("Core temperature t0:", DICTFIT["t0"])
            print("Pedestal temperature tped:", DICTFIT["tped"])
            print("Temperature at separatrix tsep:", DICTFIT["tsep"])
            print("Temperature peaking parameter alphat:", LSQFIT[0][0])
            print("Second temperature exponent betat:", LSQFIT[0][1], "\n")

            RHOARR = linspace(0, 1, 100)
            RHOPEDT = LSQFIT[0][2]
            TCORE = DICTFIT["t0"]
            TPED = DICTFIT["tped"]
            TSEP = DICTFIT["tsep"]
            ALPHAT = LSQFIT[0][0]
            BETAT = LSQFIT[0][1]

            figure()
            xlabel(r"$\rho$")
            ylabel(r"$T$")
            axvline(RHOPEDT, label=r"$\rho_{ped}$")
            axhline(TCORE, ls="--")
            axhline(TPED, ls=":")
            axhline(TSEP)
            plot(
                RHOARR,
                [
                    tprofile(x, RHOPEDT, TCORE, TPED, TSEP, ALPHAT, BETAT)
                    for x in RHOARR
                ],
            )
            plot(RHO, TEMP, "kx")


show()
