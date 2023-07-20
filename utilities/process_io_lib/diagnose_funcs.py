"""
A list of functions for the diagnose_process.py program

Author: Hanni Lux (Hanni.Lux@ccfe.ac.uk)

Notes:
19/08/2014 HL Initial version of this module.


Compatible with PROCESS version 319
"""


from process.io.mfile import MFile
from pylab import figure, plot, xticks, subplots_adjust, grid


def plot_normalised_ixc(mfilename="MFILE.DAT"):

    """
    Plots the normalised values of the iteration variables of a given MFILE
    """

    m_file = MFile(mfilename)

    list_nitvar = []
    list_labels = []

    nvar = int(m_file.data["nvar"].get_scan(-1))
    for i in range(1, nvar + 1):

        nitvar = m_file.data["nitvar{:03}".format(i)].get_scan(-1)
        list_nitvar += [nitvar]

        label = m_file.data["nitvar{:03}".format(i)].var_description

        list_labels += [label.replace("_(range_normalised)", "")]

    figure()
    subplots_adjust(bottom=0.2)
    plot(list_nitvar, "ks")
    xticks(range(len(list_nitvar)), list_labels, rotation="vertical")
    grid(True)


def plot_normalised_icc_res(mfilename="MFILE.DAT"):

    """
    Plots the normalised values of the costraint residuals of a given MFILE
    """

    m_file = MFile(mfilename)

    list_normres = []
    list_labels = []

    neqns = int(m_file.data["neqns+nineqns"].get_scan(-1))
    for i in range(1, neqns + 1):

        normres = m_file.data["normres{:03}".format(i)].get_scan(-1)
        list_normres += [normres]

        label = m_file.data["normres{:03}".format(i)].var_description

        list_labels += [label.replace("_normalised_residue", "")]

    figure()
    subplots_adjust(bottom=0.55)
    plot(list_normres, "ks")
    xticks(range(len(list_normres)), list_labels, rotation="vertical")
    grid(True)
