#!/usr/bin/env python

import sys
import argparse
import numpy as np
import matplotlib

if sys.platform == "darwin":
    matplotlib.use("WxAgg")
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
page1 = plt.figure(figsize=(6, 9))
page1.subplots_adjust(hspace=0.20, wspace=0.0)

# Second page plots
page2 = plt.figure(figsize=(6, 9))
page2.subplots_adjust(hspace=0.20, wspace=0.0)

# find max x||b
x_max = float(per_column[1][-1])

# Create numpy array of connection length
xpar = np.array([float(x) for x in per_column[1]])
# xpar[0] = 99999

# steps
steps = [int(x) for x in per_column[0]]

density = [float(x) for x in per_column[3]]
target_temperature = [float(x) for x in per_column[2]]
plasma_thermal_pressure = [float(x) for x in per_column[4]]
plasma_pressure = [float(x) for x in per_column[5]]
mach = [float(x) for x in per_column[7]]
total_neutral_density = [float(x) for x in per_column[8]]

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

# Row 1
p1r1 = page1.add_subplot(411)
p1r1.loglog(xpar, hrad_mw, label="H rad", ls="dashed")
p1r1.loglog(xpar, cx_mw, label="CX", ls="dotted")
p1r1.loglog(xpar, im_mw, label="Imp rad", ls="dashdot")
p1r1.loglog(xpar, ion_mw, label="Ionisation", ls="solid")
p1r1.set_xlim(xmin=0.0002)
p1r1.set_xlim(xmax=x_max)
p1r1.set_ylim(ymin=1.0)
# p1r1.xaxis.set_major_formatter(nullfmt)
p1r1.set_ylabel("power dens. (MWm$^{-3}$)")
p1r1.legend(loc=1, prop={"size": 10})

# # Row 2,
p1r2 = page1.add_subplot(412)
p1r2.semilogx(
    xpar, plasma_pressure, label="Plasma pressure including kinetic energy term"
)
p1r2.semilogx(xpar, plasma_thermal_pressure, label="Plasma thermal pressure")
p1r2.set_xlim(xmin=0.0002)
p1r2.set_xlim(xmax=x_max)
p1r2.set_ylim(ymin=0)
p1r2.set_ylabel("pressure (Pa)")
p1r2.legend(loc=4, prop={"size": 10})
p1r2.tick_params(axis="y", labelsize="9")
# p1r2.xaxis.set_major_formatter(nullfmt)

# Row 3,
p1r3 = page1.add_subplot(413)
p1r3.loglog(xpar, target_temperature, label="Temperature $T_e$")
p1r3.set_xlim(xmin=0.0002)
p1r3.set_xlim(xmax=x_max)
p1r3.set_ylim(ymin=1.0)
# p1r3.xaxis.set_major_formatter(nullfmt)
p1r3.set_ylabel("(eV)")
p1r3.legend(loc=4, prop={"size": 10})

# Row 4,
p1r4 = page1.add_subplot(414)
p1r4.loglog(xpar, density, label="$n_e/10^{20}m^{-3}$")
p1r4.loglog(xpar, total_neutral_density, label="$n_0/10^{20}m^{-3}$")
p1r4.loglog(xpar, n01, label="$n_{01}/10^{20}m^{-3}$", ls="dashed")
p1r4.loglog(xpar, n02, label="$n_{02}/10^{20}m^{-3}$", ls="dashed")
p1r4.loglog(xpar, mach, label="Mach")
p1r4.set_xlim(xmin=0.0002)
p1r4.set_xlim(xmax=x_max)
p1r4.set_ylim(ymin=0.01)
p1r4.tick_params(axis="x", labelsize="11")
p1r4.set_xlabel(r"Connection length from target $x_{\parallel}$ (m)", fontsize=10)
p1r4.legend(loc=4, prop={"size": 10})

# Second page plots------------------------------------
# Row 1
# Left
p2r1 = page2.add_subplot(411)
p2r1.semilogx(
    xpar,
    power_loss_integral,
    label="Integrated power emission from" " SOL\n=radiation + charge exchange",
)
p2r1.set_xlim(xmin=0.0002)
p2r1.set_xlim(xmax=x_max)
ymax = power_loss_integral[-1]
ymax = round(ymax / 50 + 0.5) * 50
# p2r1.set_ylim([0, ymax])
p2r1.set_ylabel("(MW)")
# p2r1.xaxis.set_major_formatter(nullfmt)
p2r1.legend(loc=2, prop={"size": 10})

# Row 2#
p2r2 = page2.add_subplot(412)
p2r2.set_xlim(xmin=0.0002)
p2r2.set_xlim(xmax=x_max)
p2r2.set_ylim([0.1, 1000])
p2r2.set_ylabel("power dens. (MWm$^{-3}$)")
if max(He_mw) > 0.001:
    p2r2.loglog(xpar, He_mw, label="He", ls="solid")
if max(Be_mw) > 0.001:
    p2r2.loglog(xpar, Be_mw, label="Be", ls="dashed")
if max(C_mw) > 0.001:
    p2r2.loglog(xpar, C_mw, label="C", ls="dashdot")
if max(N_mw) > 0.001:
    p2r2.loglog(xpar, N_mw, label="N", ls="dotted")
if max(O_mw) > 0.001:
    p2r2.loglog(xpar, O_mw, label="O", ls="solid")
if max(Ne_mw) > 0.001:
    p2r2.loglog(xpar, Ne_mw, label="Ne", ls="dashed")
if max(Si_mw) > 0.001:
    p2r2.loglog(xpar, Si_mw, label="Si", ls="dashdot")
if max(Ar_mw) > 0.001:
    p2r2.loglog(xpar, Ar_mw, label="Ar", ls="dotted")
if max(Fe_mw) > 0.001:
    p2r2.loglog(xpar, Fe_mw, label="Fe", ls="solid")
if max(Ni_mw) > 0.001:
    p2r2.loglog(xpar, Ni_mw, label="Ni", ls="dashed")
if max(Kr_mw) > 0.001:
    p2r2.loglog(xpar, Kr_mw, label="Kr", ls="solid")
if max(Xe_mw) > 0.001:
    p2r2.loglog(xpar, Xe_mw, label="Xe", ls="dashed")
if max(W_mw) > 0.001:
    p2r2.loglog(xpar, W_mw, label="W", ls="dashdot")
p2r2.legend(loc=1, prop={"size": 10})
p2r2.plot((x_max, x_max), (1, 10000), ls="dashed", color="black")
# p2r2.xaxis.set_major_formatter(nullfmt)

# Row 3
p2r3 = page2.add_subplot(413)
p2r3.semilogx(xpar, nv, label="Plasma flux [$10^{24}m^{-2}s^{-1}$]")
p2r3.set_xlim(xmin=0.0002)
p2r3.set_xlim(xmax=x_max)
# p2r3.xaxis.set_major_formatter(nullfmt)
p2r3.legend(loc=1, prop={"size": 10})

# Row 4
p2r4 = page2.add_subplot(414)
p2r4.semilogx(xpar, v, label="Plasma speed [ms$^{-1}$]")
p2r4.set_xlim(xmin=0.0002)
p2r4.set_xlim(xmax=x_max)
p2r4.tick_params(axis="y", labelsize="9")
p2r4.tick_params(axis="x", labelsize="11")
p2r4.set_xlabel(r"Connection length from target $x_{\parallel}$ (m)", fontsize=10)
p2r4.legend(loc=1, prop={"size": 10})

# Save as a single two-page file
with bpdf.PdfPages("1D profiles " + filename + ".pdf") as pdf:
    pdf.savefig(page1)
    pdf.savefig(page2)

page1.savefig("1D profiles " + filename + "p1.png")
page2.savefig("1D profiles " + filename + "p2.png")

page1.savefig("1D profiles " + filename + "p1.eps")
page2.savefig("1D profiles " + filename + "p2.eps")
# with bpdf.PdfPages("1D profiles " + filename + "p1.pdf") as pdf:
#         pdf.savefig(page1)
# with bpdf.PdfPages("1D profiles " + filename + "p2.pdf") as pdf:
#         pdf.savefig(page2)

plt.show(page1)
plt.show(page2)
