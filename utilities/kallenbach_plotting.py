#!/usr/bin/env python

import argparse
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf as bpdf
from matplotlib.ticker import NullFormatter

nullfmt = NullFormatter()  # no labels

parser = argparse.ArgumentParser(
    description="Plot plasma and neutral profiles \
                      in SOL for Kallenbach 1D divertor model"
)

parser.add_argument("-f", metavar="f", type=str, help="File to read")

args = parser.parse_args()

# If user has specified a filename,
if args.f:
    filename = args.f
else:
    filename = "output_divertor.txt"

f = open(filename, "r")
lines = f.readlines()
n = len(lines[8:])

# get headings
headings = list(item.strip("\n") for item in lines[7].split(" ") if item != "")
m = len(headings)
# print(headings)

per_row = []
for line in lines[8:]:
    per_row.append([it.strip("\n") for it in line.split(" ") if it != ""])

per_column = list(zip(*per_row))

# Plotting
# 0 - steps
# 1 - x//b
# 2 - te
# 3 - ne
# 4 - Ptherm
# 5 - Ptotal
# 6 - v
# 7 - mach
# 8 - n0
# 9 - power
# 10 - perparea
# 11 - qtot
# 12 - qconv
# 13 - qcond
# 14 - CX
# 15 - Ion
# 16 - Hrad
# 17 - imrad
# 18 - y7
# 19 - y8
# 20 - y9
# 21 - y10
# 22 - He
# 23 - Be
# 24 - C
# 25 - N
# 26 - O
# 27 - Ne
# 28 - Si
# 29 - Ar
# 30 - Fe
# 31 - Ni
# 32 - Kr
# 33 - Xe
# 34 - W
# 35 - n01/1e20
# 36 - n02/1e20
# 37 - nv24
# 38 - v

# First page
page1 = plt.figure(figsize=(10, 9))
page1.subplots_adjust(hspace=0.05, wspace=0.0)

# Second page plots
page2 = plt.figure(figsize=(10, 9))
page2.subplots_adjust(hspace=0.05, wspace=0.0)

# find max x||b
x_max = float(per_column[1][-1])

# Create numpy array of connection length
xpar = np.array([float(x) for x in per_column[1]])
xpar[0] = 99999


# convert to MW
cx_mw = [float(x) / 1e6 for x in per_column[14]]
ion_mw = [float(x) / 1e6 for x in per_column[15]]
hrad_mw = [float(x) / 1e6 for x in per_column[16]]
im_mw = [float(x) / 1e6 for x in per_column[17]]
# Power loss integrals are already in MW
y7_mw = np.array([float(x) for x in per_column[18]])
y8_mw = np.array([float(x) for x in per_column[19]])
y9_mw = np.array([float(x) for x in per_column[20]])
y10_mw = np.array([float(x) for x in per_column[21]])

He_mw = [float(x) / 1e6 for x in per_column[22]]
Be_mw = [float(x) / 1e6 for x in per_column[23]]
C_mw = [float(x) / 1e6 for x in per_column[24]]
N_mw = [float(x) / 1e6 for x in per_column[25]]
O_mw = [float(x) / 1e6 for x in per_column[26]]
Ne_mw = [float(x) / 1e6 for x in per_column[27]]
Si_mw = [float(x) / 1e6 for x in per_column[28]]
Ar_mw = [float(x) / 1e6 for x in per_column[29]]
Fe_mw = [float(x) / 1e6 for x in per_column[30]]
Ni_mw = [float(x) / 1e6 for x in per_column[31]]
Kr_mw = [float(x) / 1e6 for x in per_column[32]]
Xe_mw = [float(x) / 1e6 for x in per_column[33]]
W_mw = [float(x) / 1e6 for x in per_column[34]]
n01 = [float(x) for x in per_column[35]]  # /1e20 m-3
n02 = [float(x) for x in per_column[36]]  # /1e20 m-3
# Particle flux and velocity are both negative in the code
nv = [-float(x) for x in per_column[37]]  # /1e24 sm-2
v = [-float(x) for x in per_column[38]]  # ms-1

# Total power emitted by the SOL does not include the ionisation loss
power_loss_integral = y7_mw + y8_mw + y9_mw
spherical_power_load = power_loss_integral / (4 * 3.142 * (xpar * 0.5) ** 2)

# Left hand plots (odd numbers)
# Row 1
p1r1left = page1.add_subplot(421)
p1r1left.semilogy(per_column[1], hrad_mw, label="H rad", ls="dashed")
p1r1left.semilogy(per_column[1], cx_mw, label="CX", ls="dotted")
p1r1left.semilogy(per_column[1], im_mw, label="Imp rad", ls="dashdot")
p1r1left.semilogy(per_column[1], ion_mw, label="Ionisation", ls="solid")
p1r1left.set_xlim([0.0, 0.015])
# p1r1left.set_ylim(ymin=1)
p1r1left.xaxis.set_major_formatter(nullfmt)
p1r1left.set_ylabel("power dens. (MWm$^{-3}$)")


# Row 2, Left
p1r2left = page1.add_subplot(423)
p1r2left.plot(
    per_column[1], per_column[5], label="Plasma pressure including kinetic energy term"
)
p1r2left.plot(per_column[1], per_column[4], label="Plasma thermal pressure")
p1r2left.set_xlim([0.0, 0.015])
p1r2left.set_ylim(ymin=0)
# p1r2left.set_ylim(ymax=2700)
p1r2left.set_ylabel("pressure (Pa)")
# p1r2left.legend(loc=4, prop={'size': 12})
p1r2left.tick_params(axis="y", labelsize="9")
p1r2left.xaxis.set_major_formatter(nullfmt)

# Row 3, Left
p1r3left = page1.add_subplot(425)
p1r3left.semilogy(per_column[1], per_column[2], label="$T_e$")
p1r3left.set_xlim([0.0, 0.015])
p1r3left.set_ylim(ymin=1.0)
p1r3left.xaxis.set_major_formatter(nullfmt)
p1r3left.set_ylabel("(eV)")


# Row 4, Left
p1r4left = page1.add_subplot(427)
p1r4left.semilogy(per_column[1], per_column[3], label="$n_e/10^{20}m^{-3}$")
p1r4left.semilogy(per_column[1], per_column[8], label="$n_0/10^{20}m^{-3}$")
p1r4left.semilogy(per_column[1], n01, label="$n_01/10^{20}m^{-3}$", ls="dashed")
p1r4left.semilogy(per_column[1], n02, label="$n_02/10^{20}m^{-3}$", ls="dashed")
p1r4left.semilogy(per_column[1], per_column[7], label="mach")
# p1r4left.set_ylim([0.1, 100])
p1r4left.set_xlim([0.0, 0.015])
p1r4left.set_ylim(ymin=0.01)
p1r4left.set_xlabel(r"$x\parallel B$ (m)", fontsize=14)


# Row 1, Right
p1r1right = page1.add_subplot(422)
p1r1right.loglog(per_column[1], hrad_mw, label="H rad", ls="dashed")
p1r1right.loglog(per_column[1], cx_mw, label="CX", ls="dotted")
p1r1right.loglog(per_column[1], im_mw, label="Imp rad", ls="dashdot")
p1r1right.loglog(per_column[1], ion_mw, label="Ionisation", ls="solid")
p1r1right.set_xlim([0.015, x_max])
# p1r1right.set_ylim(ymin=1)
p1r1right.legend(loc=1, prop={"size": 10})
# p1r1right.plot((x_max, x_max), (0.001, 10000), ls='dashed', color="black")
# no labels
p1r1right.xaxis.set_major_formatter(nullfmt)
p1r1right.yaxis.set_major_formatter(nullfmt)
p1r1right.set_ylim(p1r1left.get_ylim())

# Row 2, Right
p1r2right = page1.add_subplot(424)
p1r2right.semilogx(
    per_column[1], per_column[5], label="Plasma pressure including kinetic energy term"
)
p1r2right.semilogx(per_column[1], per_column[4], label="Plasma thermal pressure")
p1r2right.set_xlim([0.015, x_max])
p1r2right.set_ylim(ymin=0)
p1r2right.legend(loc=4, prop={"size": 10})
p1r2right.tick_params(axis="y", labelsize="9")
# no labels
p1r2right.xaxis.set_major_formatter(nullfmt)
p1r2right.yaxis.set_major_formatter(nullfmt)
p1r2right.set_ylim(p1r2left.get_ylim())

# Row 3, Right
p1r3right = page1.add_subplot(426)
p1r3right.loglog(per_column[1], per_column[2], label="$T_e$")
p1r3right.set_xlim([0.015, x_max])
# p1r3right.set_ylim([0.0, 350.0])
p1r3right.set_ylim(ymin=1.0)
p1r3right.legend(loc=4, prop={"size": 10})
# no labels
p1r3right.xaxis.set_major_formatter(nullfmt)
p1r3right.yaxis.set_major_formatter(nullfmt)
p1r3right.set_ylim(p1r3left.get_ylim())

# Row 4, Right
p1r4right = page1.add_subplot(428)
p1r4right.loglog(per_column[1], per_column[3], label="$n_e/10^{20}m^{-3}$")
p1r4right.loglog(per_column[1], per_column[8], label="$n_0/10^{20}m^{-3}$")
p1r4right.loglog(per_column[1], n01, label="$n_{01}/10^{20}m^{-3}$", ls="dashed")
p1r4right.loglog(per_column[1], n02, label="$n_{02}/10^{20}m^{-3}$", ls="dashed")
p1r4right.loglog(per_column[1], per_column[7], label="mach")
# p1r4right.set_ylim([0.1, 10])
p1r4right.set_xlim([0.015, x_max])
# p1r4right.set_ylim(ymin=0.1)
p1r4right.set_xlabel(r"$x\parallel B$ (m)", fontsize=14)
p1r4right.legend(loc=1, prop={"size": 10})
# no labels
p1r4right.yaxis.set_major_formatter(nullfmt)
p1r4right.set_ylim(p1r4left.get_ylim())

# Second page plots
# Row 1
# Left
p2r1left = page2.add_subplot(421)
p2r1left.plot(
    per_column[1], power_loss_integral, label="Integrated power emission from SOL"
)
p2r1left.set_xlim([0.0, 0.015])
ymax = power_loss_integral[-1]
ymax = round(ymax / 50 + 0.5) * 50
# p2r1left.set_ylim([0, ymax])
p2r1left.set_ylabel("(MW)")
p2r1left.xaxis.set_major_formatter(nullfmt)

# Right
p2r1right = page2.add_subplot(422)
p2r1right.semilogx(
    per_column[1],
    power_loss_integral,
    label="Integrated power emission from SOL\n=radiation + charge exchange",
)

p2r1right.set_xlim([0.015, x_max])
# p2r1right.set_ylim([0, ymax])
# p2r1right.set_xlabel("$x\parallel B$ (m)", fontsize=14)
p2r1right.legend(loc=4, prop={"size": 10})
# no labels
p2r1right.xaxis.set_major_formatter(nullfmt)
p2r1right.yaxis.set_major_formatter(nullfmt)
p2r1right.set_ylim(p2r1left.get_ylim())

# Row 2
# Left
p2r2left = page2.add_subplot(423)
p2r2left.set_xlim([0.0, 0.015])
p2r2left.set_ylim([0.1, 1000])
# p2r2left.set_xlabel("$x\parallel B$ (m)", fontsize=14)
p2r2left.set_ylabel("power dens. (MWm$^{-3}$)")
if max(He_mw) > 0.001:
    p2r2left.semilogy(per_column[1], He_mw, label="He", ls="solid")
if max(Be_mw) > 0.001:
    p2r2left.semilogy(per_column[1], Be_mw, label="Be", ls="dashed")
if max(C_mw) > 0.001:
    p2r2left.semilogy(per_column[1], C_mw, label="C", ls="dashdot")
if max(N_mw) > 0.001:
    p2r2left.semilogy(per_column[1], N_mw, label="N", ls="dotted")
if max(O_mw) > 0.001:
    p2r2left.semilogy(per_column[1], O_mw, label="O", ls="solid")
if max(Ne_mw) > 0.001:
    p2r2left.semilogy(per_column[1], Ne_mw, label="Ne", ls="dashed")
if max(Si_mw) > 0.001:
    p2r2left.semilogy(per_column[1], Si_mw, label="Si", ls="dashdot")
if max(Ar_mw) > 0.001:
    p2r2left.semilogy(per_column[1], Ar_mw, label="Ar", ls="dotted")
if max(Fe_mw) > 0.001:
    p2r2left.semilogy(per_column[1], Fe_mw, label="Fe", ls="solid")
if max(Ni_mw) > 0.001:
    p2r2left.semilogy(per_column[1], Ni_mw, label="Ni", ls="dashed")
if max(Kr_mw) > 0.001:
    p2r2left.semilogy(per_column[1], Kr_mw, label="Kr", ls="solid")
if max(Xe_mw) > 0.001:
    p2r2left.semilogy(per_column[1], Xe_mw, label="Xe", ls="dashed")
if max(W_mw) > 0.001:
    p2r2left.semilogy(per_column[1], W_mw, label="W", ls="dashdot")
# p2r2left.legend(loc=1, prop={'size': 8})
p2r2left.plot((x_max, x_max), (1, 10000), ls="dashed", color="black")
p2r2left.xaxis.set_major_formatter(nullfmt)


# Row 2 Right
p2r2right = page2.add_subplot(424)
p2r2right.set_xlim([0.015, x_max])
p2r2right.set_ylim([0.1, 1000])
# p2r2right.set_xlabel("$x\parallel B$ (m)", fontsize=14)
# p2r2right.set_ylabel("power dens. (MWm$^{-3}$)")
if max(He_mw) > 0.001:
    p2r2right.loglog(per_column[1], He_mw, label="He", ls="solid")
if max(Be_mw) > 0.001:
    p2r2right.loglog(per_column[1], Be_mw, label="Be", ls="dashed")
if max(C_mw) > 0.001:
    p2r2right.loglog(per_column[1], C_mw, label="C", ls="dashdot")
if max(N_mw) > 0.001:
    p2r2right.loglog(per_column[1], N_mw, label="N", ls="dotted")
if max(O_mw) > 0.001:
    p2r2right.loglog(per_column[1], O_mw, label="O", ls="solid")
if max(Ne_mw) > 0.001:
    p2r2right.loglog(per_column[1], Ne_mw, label="Ne", ls="dashed")
if max(Si_mw) > 0.001:
    p2r2right.loglog(per_column[1], Si_mw, label="Si", ls="dashdot")
if max(Ar_mw) > 0.001:
    p2r2right.loglog(per_column[1], Ar_mw, label="Ar", ls="dotted")
if max(Fe_mw) > 0.001:
    p2r2right.loglog(per_column[1], Fe_mw, label="Fe", ls="solid")
if max(Ni_mw) > 0.001:
    p2r2right.loglog(per_column[1], Ni_mw, label="Ni", ls="dashed")
if max(Kr_mw) > 0.001:
    p2r2right.loglog(per_column[1], Kr_mw, label="Kr", ls="solid")
if max(Xe_mw) > 0.001:
    p2r2right.loglog(per_column[1], Xe_mw, label="Xe", ls="dashed")
if max(W_mw) > 0.001:
    p2r2right.loglog(per_column[1], W_mw, label="W", ls="dashdot")
p2r2right.legend(loc=1, prop={"size": 10})
# no labels
p2r2right.xaxis.set_major_formatter(nullfmt)
p2r2right.yaxis.set_major_formatter(nullfmt)
p2r2right.set_ylim(p2r2left.get_ylim())

# Row 3 Left
p2r3left = page2.add_subplot(425)
p2r3left.plot(per_column[1], nv, label="Plasma flux [$10^{24}m^{-2}s^{-1}$]")
p2r3left.set_xlim([0.0, 0.015])
p2r3left.set_ylim(ymax=9.9)
ylim = p2r3left.get_ylim()

p2r3left.xaxis.set_major_formatter(nullfmt)
# ymax = v[-1]
# ymax = round(ymax/50 + 0.5) * 50
# p2r3left.set_ylim([0, ymax])
# p2r3left.legend(loc=1, prop={'size': 8})

# Row 3 Right
p2r3right = page2.add_subplot(426)
p2r3right.semilogx(per_column[1], nv, label="Plasma flux [$10^{24}m^{-2}s^{-1}$]")
p2r3right.set_xlim([0.015, x_max])
p2r3right.set_ylim([0, ylim[1]])
p2r3right.legend(loc=1, prop={"size": 10})
# no labels
p2r3right.xaxis.set_major_formatter(nullfmt)
p2r3right.yaxis.set_major_formatter(nullfmt)
p2r3right.set_ylim(p2r3left.get_ylim())

# Row 4 Left
p2r4left = page2.add_subplot(427)
p2r4left.plot(per_column[1], v, label="Plasma speed [ms$^{-1}$]")
p2r4left.set_xlim([0.0, 0.015])
ylim = p2r4left.get_ylim()
p2r4left.tick_params(axis="y", labelsize="9")
p2r4left.set_xlabel(r"$x\parallel B$ (m)", fontsize=14)

# Row 4 Right
p2r4right = page2.add_subplot(428)
p2r4right.semilogx(per_column[1], v, label="Plasma speed [ms$^{-1}$]")
p2r4right.set_xlim([0.015, x_max])
p2r4right.set_xlabel(r"$x\parallel B$ (m)", fontsize=14)
p2r4right.set_ylim([0, ylim[1]])
p2r4right.legend(loc=1, prop={"size": 10})
p2r4right.tick_params(axis="y", labelsize="9")
# no labels
p2r4right.yaxis.set_major_formatter(nullfmt)
p2r4right.set_ylim(p2r4left.get_ylim())

# This should be called after all axes have been added
# Unfortunately it seems to override the hspace and wspace parameters.
# page1.tight_layout()
# page2.tight_layout()

# Save as a single two-page file
with bpdf.PdfPages("1D profiles " + filename + ".pdf") as pdf:
    pdf.savefig(page1)
    pdf.savefig(page2)

with bpdf.PdfPages("1D profiles " + filename + "p1.pdf") as pdf:
    pdf.savefig(page1)
with bpdf.PdfPages("1D profiles " + filename + "p2.pdf") as pdf:
    pdf.savefig(page2)

plt.show(page1)
plt.show(page2)
