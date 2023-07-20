#!/usr/bin/env python
"""

  Plots profiles from PROCESS for a number of MFILEs

  James Morris
  15/02/2019
  CCFE

"""

import argparse
import process.io.mfile as mf
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf as bpdf
import numpy as np

matplotlib.rcParams["figure.max_open_warning"] = 40


def plot_nprofile(prof, mfdat, scan):
    """Function to plot density profile
    Arguments:
      prof --> axis object to add plot to
      mfdat --> MFILE data object
      scan --> scan number to plot
    """
    xmin = 0
    xmax = 1
    ymin = 0
    ymax = 20
    prof.set_ylim([ymin, ymax])
    prof.set_xlim([xmin, xmax])
    prof.set_autoscaley_on(False)
    prof.set_xlabel("r/a")
    prof.set_ylabel("ne / 1e19 m-3")
    prof.set_title("Density profile")

    ipedestal = mfdat.data["ipedestal"].get_scan(scan)
    neped = mfdat.data["neped"].get_scan(scan)
    nesep = mfdat.data["nesep"].get_scan(scan)
    rhopedn = mfdat.data["rhopedn"].get_scan(scan)
    alphan = mfdat.data["alphan"].get_scan(scan)
    ne0 = mfdat.data["ne0"].get_scan(scan)

    if ipedestal == 1:
        rhocore1 = np.linspace(0, 0.95 * rhopedn)
        rhocore2 = np.linspace(0.95 * rhopedn, rhopedn)
        rhocore = np.append(rhocore1, rhocore2)
        ncore = neped + (ne0 - neped) * (1 - rhocore**2 / rhopedn**2) ** alphan

        rhosep = np.linspace(rhopedn, 1)
        nsep = nesep + (neped - nesep) * (1 - rhosep) / (1 - min(0.9999, rhopedn))

        rho = np.append(rhocore, rhosep)
        ne = np.append(ncore, nsep)
    else:
        rho1 = np.linspace(0, 0.95)
        rho2 = np.linspace(0.95, 1)
        rho = np.append(rho1, rho2)
        ne = ne0 * (1 - rho**2) ** alphan
    ne = ne / 1e19
    prof.plot(rho, ne, label=mfdat.filename)


def plot_tprofile(prof, mfdat, scan):
    """Function to plot temperature profile
    Arguments:
      prof --> axis object to add plot to
      mfdat --> MFILE data object
      scan --> scan number to plot
    """
    xmin = 0
    xmax = 1
    ymin = 0
    ymax = 40
    prof.set_ylim([ymin, ymax])
    prof.set_xlim([xmin, xmax])
    prof.set_autoscaley_on(False)
    prof.set_xlabel("r/a")
    prof.set_ylabel("Te / keV")
    prof.set_title("Temperature profile")

    ipedestal = mfdat.data["ipedestal"].get_scan(scan)
    rhopedt = mfdat.data["rhopedt"].get_scan(scan)
    tbeta = mfdat.data["tbeta"].get_scan(scan)
    teped = mfdat.data["teped"].get_scan(scan)
    tesep = mfdat.data["tesep"].get_scan(scan)
    alphat = mfdat.data["alphat"].get_scan(scan)
    te0 = mfdat.data["te0"].get_scan(scan)

    if ipedestal == 1:
        rhocore1 = np.linspace(0, 0.9 * rhopedt)
        rhocore2 = np.linspace(0.9 * rhopedt, rhopedt)
        rhocore = np.append(rhocore1, rhocore2)
        tcore = teped + (te0 - teped) * (1 - (rhocore / rhopedt) ** tbeta) ** alphat

        rhosep = np.linspace(rhopedt, 1)
        tsep = tesep + (teped - tesep) * (1 - rhosep) / (1 - min(0.9999, rhopedt))

        rho = np.append(rhocore, rhosep)
        te = np.append(tcore, tsep)
    else:
        rho = np.linspace(0, 1)
        te = te0 * (1 - rho**2) ** alphat
    prof.plot(rho, te, label=mfdat.filename)


if __name__ == "__main__":

    # Setup command line arguments
    parser = argparse.ArgumentParser(
        description="Produces plots of density and temperature profiles, using MFILE(s).  "
        "For info contact james.morris2@ukaea.uk"
    )

    parser.add_argument(
        "-f",
        metavar="files",
        type=str,
        nargs="+",
        default="MFILE.DAT",
        help='specify input/output file paths as a list "-f ../file1 ../file2"',
    )

    parser.add_argument(
        "-s", "--save", help="save plot as well as showing figure", action="store_true"
    )

    parser.add_argument("-n", type=int, help="Which scan to plot?")

    parser.add_argument(
        "-o",
        metavar="output_name",
        type=str,
        default="profiles.pdf",
        help="name of output pdf",
    )

    args = parser.parse_args()

    if args.n:
        scan = args.n
    else:
        scan = -1

    # Setup plot area
    fig = plt.figure(figsize=(12, 10), dpi=80)

    # Density plot
    plot_1 = fig.add_subplot(211, aspect=1 / 20)

    # Temperature plot
    plot_2 = fig.add_subplot(212, aspect=1 / 40)

    # read MFILE
    if args.f != "MFILE.DAT":
        for item in args.f:
            m_file = mf.MFile(item)
            plot_nprofile(plot_1, m_file, scan)
            plot_tprofile(plot_2, m_file, scan)

    else:
        m_file = mf.MFile("MFILE.DAT")
        m_file = mf.MFile(item)
        plot_nprofile(plot_1, m_file, scan)
        plot_tprofile(plot_2, m_file, scan)

    plot_1.legend(loc="upper left", bbox_to_anchor=(1, 1))
    plot_2.legend(loc="upper left", bbox_to_anchor=(1, 1))

    if args.save:
        with bpdf.PdfPages(args.o) as pdf:
            pdf.savefig(fig)

    plt.show(fig)
